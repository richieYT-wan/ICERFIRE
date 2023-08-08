from datetime import datetime as dt
import os, sys
import pandas as pd
import argparse

module_path = os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
    sys.path.append(module_path)
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/'
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

# Get the absolute path of the "src" directory and HARDCODEDDLY ADDING IT
src_path = '/home/local/tools/src/ICERFIRE-1.0/src/'
# Add the "src" directory to the Python module search path
sys.path.append(src_path)
# print(module_path, sys.path)
from train_eval import evaluate_trained_models
from mutation_tools import pipeline_mutation_scores
from utils import str2bool, get_random_id, get_datetime_string, mkdirs, pkl_load


def get_rank(pred, hp):
    return (hp[0] > pred).mean() * 100


def args_parser():
    parser = argparse.ArgumentParser(
        'Runs the ICERFIRE model on preprocessed data, assuming data has been processed to return' \
        'NetMHCpan ranks, ICOREs, self-similarity score, and PepX expression scores (Optional).')
    parser.add_argument('-j', '--jobid', dest='jobid', type=str, help='Job ID from the server')
    parser.add_argument('-f', '--infile', dest='infile', type=str, required=True,
                        help='Full path to the file containing the icores/selfsimilarity scores')
    parser.add_argument('-pf', '--pepxpath', dest='pepxpath', type=str, required=False, default=None,
                        help='Full path to the file containing the PepX query of the test file')
    parser.add_argument('-ae', '--add_expression', dest='add_expression', type=str2bool,
                        required=False, default=True,
                        help='Whether to use the model that includes expression as a feature')
    parser.add_argument('-ue', '--user_expression', dest='user_expression', type=str2bool,
                        default=False, help='Whether the user provides their own expression values')
    parser.add_argument('-o', '--outdir', dest='outdir', type=str, required=False, default='../tmp/',
                        help='Output directory')
    return parser.parse_args()


def main():
    args = vars(args_parser())

    # Get the output directory with a random ID and a tag to indicate whether we used the model with expression
    run_dt = get_datetime_string()
    run_id = get_random_id(6)
    run_tag = 'AddExpr' if args['add_expression'] else 'NoExpr'
    # basename = os.path.basename(args['infile']).split('.')[0]
    # print('basename', basename)
    # print('infile, pepxfile', args['infile'], args['pepxpath'])
    run_name = f'{run_dt}_{run_tag}_{run_id}'
    outdir = os.path.join(args['outdir'], f'{run_name}/')
    mkdirs(outdir)
    jobid = str(args['jobid'])
    # Get the directory one level above the script
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/'
    # Load appropriate model and data
    data = pd.read_csv(args['infile'], sep=' ')
    preds_100k = pd.read_csv('/tools/src/ICERFIRE-1.0/data/human_proteome/hp_preds_100k.txt', header=None)
    unpickle = pkl_load(f'{parent_dir}saved_models/ICERFIRE_Expr{args["add_expression"]}.pkl')
    models, kwargs, ics = unpickle['model'], unpickle['kwargs'], unpickle['ics']

    if args['add_expression'] and args['user_expression'] and 'total_gene_tpm' not in data.columns:
        args['add_expression'] = False
        unpickle = pkl_load(f'{parent_dir}saved_models/ICERFIRE_Expr{False}.pkl')
        models, kwargs, ics = unpickle['model'], unpickle['kwargs'], unpickle['ics']
        print('User-provided expression was selected but no TPM values found in the data. Continuing with a model without expression')

    if args['add_expression'] and os.path.exists(args['pepxpath']) and args['pepxpath'] != "None":
        # TODO : DEAL WITH case where PepX is not used and maybe expression is still enabled (and provided)
        pepx = pd.read_csv(args['pepxpath'])
        data = pd.merge(data, pepx.rename(columns={'peptide': 'icore_wt_aligned'}), how='left',
                        left_on='icore_wt_aligned', right_on='icore_wt_aligned')
        data["TPMFilledWithMedian"] = data["total_gene_tpm"].isna()
        median_value = data["total_gene_tpm"].median(skipna=True)
        data["total_gene_tpm"].fillna(median_value, inplace=True)

        data.fillna(data.median(skipna=True, numeric_only=True), inplace=True)

    data = pipeline_mutation_scores(data, 'icore_mut', 'icore_wt_aligned', ics,
                                    threshold=kwargs['threshold'], prefix='icore_')
    data['seq_id'] = [f'seq_{i}' for i in range(1, len(data) + 1)]

    predictions, test_results = evaluate_trained_models(data, models, ics, encoding_kwargs=kwargs, test_mode=True,
                                                        n_jobs=8)
    # Saving results as CSV table
    predictions.sort_values('Peptide', ascending=True, inplace=True)
    predictions.rename(columns={'mean_pred': 'prediction'}, inplace=True)
    predictions.reset_index(drop=True, inplace=True)
    predictions['%Rank'] = predictions['prediction'].apply(get_rank, hp=preds_100k)

    cols_to_save = ['Peptide', 'wild_type', 'HLA', 'Pep', 'Core', 'icore_start_pos', 'icore_mut', 'icore_wt_aligned', 'EL_rank_mut',
                    'EL_rank_wt_aligned']
    cols_to_save = cols_to_save + kwargs['mut_col'] + ['prediction', '%Rank']
    predictions.to_csv(f'{outdir}ICERFIRE_predictions.csv',
                       columns=cols_to_save, index=False)
    if test_results is not None:
        pd.DataFrame(test_results).rename(columns={k: v for k, v in zip(range(len(test_results.keys())),
                                                                        [f'fold_{x}' for x in
                                                                         range(1, len(test_results.keys()))])}) \
            .to_csv(f'{outdir}ICERFIRE_metrics_per_fold.csv', index=False)

    # Cleaning input/temporary files and returning the final saved filename
    for f in os.listdir(args['outdir']):
        if f.endswith('.csv') or f.endswith('.txt'):
            os.remove(os.path.join(args['outdir'], f))
    return predictions, run_name, jobid


if __name__ == '__main__':
    predictions, run_name, jobid = main()
    print('Click ' + '<a href="https://services.healthtech.dtu.dk/services/ICERFIRE-1.0/tmp/' \
          + f'{jobid}/{run_name}/' \
            'ICERFIRE_predictions.csv" target="_blank">here</a>' + ' to download the results in .csv format.')
    pd.set_option('display.max_columns', 30)
    pd.set_option('display.max_rows', 99)
    print("\n \nBelow is a table preview of the ICERFIRE predictions as well as the identified ICOREs. \n \n")
    print(predictions[['Peptide', 'wild_type', 'icore_mut', 'icore_wt_aligned', 'EL_rank_mut', 'EL_rank_wt_aligned',
                       'total_gene_tpm', 'prediction', '%Rank']])
