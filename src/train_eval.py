import copy
import multiprocessing

import pandas as pd
import numpy as np
from src.data_processing import get_dataset, assert_encoding_kwargs
from src.metrics import get_metrics, get_mean_roc_curve, get_predictions
import sklearn
from sklearn.model_selection import ParameterGrid
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
from functools import partial
from tqdm.auto import tqdm


# TRAIN WITH PARALLEL WRAPPER
def parallel_inner_train_wrapper(train_dataframe, x_test, base_model, ics_dict,
                                 encoding_kwargs, standardize, fold_out, fold_in):
    seed = fold_out * 10 + fold_in
    # Copy the base model, resets the seed
    model = sklearn.base.clone(base_model)
    model.set_params(random_state=seed)
    if standardize:
        model = Pipeline([('scaler', StandardScaler()), ('model', model)])

    # Here resets model weight at every fold, using the fold number (range(0, n_folds*(n_folds-1)) ) as seed
    # Query subset dataframe and get encoded data and targets
    train = train_dataframe.query('fold != @fold_out and fold != @fold_in').reset_index(drop=True)
    valid = train_dataframe.query('fold == @fold_in').reset_index(drop=True)
    # Get datasets
    x_train, y_train = get_dataset(train, ics_dict, **encoding_kwargs)
    x_valid, y_valid = get_dataset(valid, ics_dict, **encoding_kwargs)

    # Fit the model and append it to the list
    model.fit(x_train, y_train)

    # Get the prediction values on both the train and validation set
    y_train_pred, y_train_score = model.predict(x_train), model.predict_proba(x_train)[:, 1]
    y_val_pred, y_val_score = model.predict(x_valid), model.predict_proba(x_valid)[:, 1]
    # Get the metrics and save them
    try:
        train_metrics = get_metrics(y_train, y_train_score, y_train_pred)
    except:
        print(train_dataframe.head())
        raise ValueError(f'{encoding_kwargs}')
    try:
        valid_metrics = get_metrics(y_valid, y_val_score, y_val_pred)
    except:
        print(train_dataframe.head())
        raise ValueError(f'{encoding_kwargs}')
    y_pred_test = model.predict_proba(x_test)[:, 1]

    return model, train_metrics, valid_metrics, y_pred_test


def nested_kcv_train(dataframe, base_model, ics_dict, encoding_kwargs: dict = None, n_jobs: int = None):
    """
    Args:
        dataframe:
        base_model:
        ics_dict:
        encoding_kwargs:
        n_jobs (int): number of parallel processes. If None, will use len(inner_folds)

    Returns:
        models_fold
        train_results
        test_results
    """
    encoding_kwargs = assert_encoding_kwargs(encoding_kwargs, mode_eval=False)
    models_dict = {}
    test_metrics = {}
    train_metrics = {}
    folds = sorted(dataframe.fold.unique())
    std = encoding_kwargs.pop('standardize')
    for fold_out in tqdm(folds, leave=False, desc='Train:Outer fold', position=2):
        # Get test set & init models list to house all models trained in inner fold
        test = dataframe.query('fold == @fold_out').reset_index(drop=True)
        x_test, y_test = get_dataset(test, ics_dict, **encoding_kwargs)
        # For a given fold, all the models that are trained should be appended to this list
        inner_folds = sorted([f for f in folds if f != fold_out])
        # N jobs must be lower than cpu_count
        n_jobs = min(multiprocessing.cpu_count() - 1, len(inner_folds)) if n_jobs is None \
            else n_jobs if (n_jobs is not None and n_jobs <= multiprocessing.cpu_count()) \
            else multiprocessing.cpu_count() - 1
        # Create the sub-dict and put the model into the models dict
        train_wrapper_ = partial(parallel_inner_train_wrapper, train_dataframe=dataframe, x_test=x_test,
                                 base_model=base_model, ics_dict=ics_dict, encoding_kwargs=encoding_kwargs,
                                 standardize=std, fold_out=fold_out)
        output = Parallel(n_jobs=n_jobs)(
            delayed(train_wrapper_)(fold_in=fold_in) for fold_in in tqdm(inner_folds,
                                                                         desc='Inner Folds',
                                                                         leave=False, position=1))
        models_dict[fold_out] = [x[0] for x in output]
        train_tmp = [x[1] for x in output]
        valid_tmp = [x[2] for x in output]
        avg_prediction = [x[3] for x in output]
        avg_prediction = np.mean(np.stack(avg_prediction), axis=0)
        train_metrics[fold_out] = {k: {'train': v_train,
                                       'valid': v_valid} for k, v_train, v_valid in
                                   zip(inner_folds, train_tmp, valid_tmp)}
        test_metrics[fold_out] = get_metrics(y_test, avg_prediction)

    return models_dict, train_metrics, test_metrics


# EVAL WITH PARALLEL WRAPPER
def parallel_eval_wrapper(test_dataframe, models_list, ics_dict, encoding_kwargs, fold_out, test_mode=False,
                          kcv_eval=False):
    # If no train dataframe provided and test_dataframe is partitioned,
    # It will eval on each of the folds
    if kcv_eval and 'fold' in test_dataframe.columns:
        test_df = test_dataframe.query('fold==@fold_out')
    else:
        test_df = test_dataframe.copy().reset_index(drop=True)

    # TODO: HERE NEED TO FIX BEHAVIOUR WHERE TARGET LABEL IS NOT PROVIDED
    # print(5, len(test_dataframe))
    # print(6, len(test_df))
    predictions_df = get_predictions(test_df, models_list, ics_dict, encoding_kwargs, test_mode).assign(fold=fold_out)
    # print(7, len(predictions_df))
    if encoding_kwargs['target_col'] in predictions_df.columns:
        test_metrics = get_metrics(predictions_df[encoding_kwargs['target_col']].values,
                                   predictions_df['pred'].values)
    else:
        test_metrics = None
    return predictions_df, test_metrics


def evaluate_trained_models(test_dataframe, models_dict, ics_dict, encoding_kwargs: dict = None, test_mode=False,
                            kcv_eval=False, n_jobs=None):
    """

    Args:
        test_mode:
        test_dataframe:
        models_dict:
        ics_dict:
        train_metrics (dict): Should be used if standardize in encoding_kwargs is True...
        encoding_kwargs:

    Returns:
        test_results
        predictions_df
    """
    encoding_kwargs = assert_encoding_kwargs(encoding_kwargs, mode_eval=True)
    # Wrapper and parallel evaluation
    eval_wrapper_ = partial(parallel_eval_wrapper, test_dataframe=test_dataframe, ics_dict=ics_dict, kcv_eval=kcv_eval,
                            encoding_kwargs=encoding_kwargs, test_mode=test_mode)
    n_jobs = len(models_dict.keys()) if (
            n_jobs is None and len(models_dict.keys()) <= multiprocessing.cpu_count()) else n_jobs
    # TODO: HERE NEED TO FIX BEHAVIOUR WHERE TARGET LABEL IS NOT PROVIDED
    output = Parallel(n_jobs=n_jobs)(delayed(eval_wrapper_)(fold_out=fold_out, models_list=models_list) \
                                     for (fold_out, models_list) in tqdm(models_dict.items(),
                                                                         desc='Eval Folds',
                                                                         leave=False,
                                                                         position=2))
    predictions_df = [x[0] for x in output]
    # print('here', len(predictions_df), len(predictions_df[0]))
    # Here simply concatenates it to get all the predictions from the folds

    predictions_df = pd.concat(predictions_df)
    print(8, len(predictions_df))
    # Get the mean predictions
    predictions_df = predictions_df.groupby(['Peptide', 'HLA', 'seq_id']).agg(mean_pred=('pred', 'mean')).reset_index()
    if encoding_kwargs['target_col'] in predictions_df.columns:
        test_metrics = [x[1] for x in output]
        test_results = {k: v for k, v in zip(models_dict.keys(), test_metrics)}
    else:
        test_results = None
    predictions_df = pd.merge(test_dataframe, predictions_df, left_on=['Peptide', 'HLA', 'seq_id'],
                              right_on=['Peptide', 'HLA', 'seq_id']) \
        .sort_values('seq_id', ascending=True).drop(columns=['seq_id'])
    return predictions_df, test_results
