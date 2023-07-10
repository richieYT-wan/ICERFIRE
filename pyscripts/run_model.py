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

# print(module_path, sys.path)
from src.utils import pkl_load
from src.train_eval import evaluate_trained_models
from src.mutation_tools import pipeline_mutation_scores
from src.utils import str2bool, get_random_id, get_datetime_string, mkdirs


def args_parser():
    parser = argparse.ArgumentParser(
        'Runs the ICERFIRE model on preprocessed data, assuming data has been processed to return' \
        'NetMHCpan ranks, ICOREs, self-similarity score, and PepX expression scores (Optional).')
    parser.add_argument('-f', '--filepath', dest='filepath', type=str, required=True,
                        help='Full path to the file containing the icores/selfsimilarity scores')
    parser.add_argument('-pf', '--pepxpath', dest='pepxpath', type=str, required=False, default=None,
                        help='Full path to the file containing the PepX query of the test file')
    parser.add_argument('-ae', '--add_expression', dest='add_expression', type=str2bool,
                        required=False, default=True,
                        help='Whether to use the model that includes expression as a feature')
    parser.add_argument('-o', '--outdir', dest='outdir', type=str, required=False, default='../tmp/',
                        help='Output directory')
    return parser.parse_args()


def main():
    args = vars(args_parser())
    # Get the output directory with a random ID and a tag to indicate whether we used the model with expression
    run_dt = get_datetime_string()
    run_id = get_random_id(6)
    run_tag = 'AddExpr' if args['add_expression'] else 'NoExpr'
    basename = os.path.basename(args['filepath'].split('.')[0])
    run_name = f'{run_dt}_{basename}_{run_tag}_{run_id}'
    outdir = os.path.join(args['outdir'], f'{run_name}/')
    mkdirs(outdir)

    # Get the directory one level above the script
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/'
    # Load appropriate model
    # TODO: TRAIN/ADD MODEL W/O EXPRESSION
    unpickle = pkl_load(f'{parent_dir}saved_models/ICERFIRE_Expr{args["add_expression"]}.pkl')
    models, kwargs, ics = unpickle['model'], unpickle['kwargs'], unpickle['ics']
    data = pd.read_csv(args['filepath'], sep=' ')
    # print(1, len(data))
    if args['add_expression'] and os.path.exists(args['pepxpath']) and args['pepxpath']!="None":
        # TODO : DEAL WITH case where PepX is not used and maybe expression is still enabled (and provided)
        pepx = pd.read_csv(args['pepxpath'])
        data = pd.merge(data, pepx.rename(columns={'peptide': 'icore_wt_aligned'}), how='left',
                        left_on='icore_wt_aligned', right_on='icore_wt_aligned')
        data.fillna(data.median(skipna=True, numeric_only=True), inplace=True)
    # print(2, len(data))

    data = pipeline_mutation_scores(data, 'icore_mut', 'icore_wt_aligned', ics,
                                    threshold=kwargs['threshold'], prefix='icore_')
    data['seq_id'] = [f'seq_{i}' for i in range(1,len(data)+1)]
    # print(3, len(data))

    predictions, test_results = evaluate_trained_models(data, models, ics, encoding_kwargs=kwargs, test_mode=True, n_jobs=8)
    # Saving results as CSV table
    # print(4, len(predictions))
    predictions.sort_values('Peptide', ascending=False).to_csv(f'{outdir}{run_name}_predictions.csv', index=False)
    if test_results is not None:
        pd.DataFrame(test_results).rename(columns={k: v for k, v in zip(range(len(test_results.keys())),
                                                                        [f'fold_{x}' for x in range(1, len(test_results.keys()))])})\
                                  .to_csv(f'{outdir}{run_name}_metrics_per_fold.csv', index=False)


if __name__ == '__main__':
    main()
