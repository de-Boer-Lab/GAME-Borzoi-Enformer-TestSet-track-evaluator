'''Calculate and save the final evaluation metrics.'''

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented

import os
import sys
import json
import pandas as pd
import numpy as np
from datetime import datetime, timezone
from scipy.stats import pearsonr
from scipy.spatial import distance
from config import EVALUATOR_NAME, EVALUATOR_INPUT_PATH, MEASURED_DATA_PATH
import pyBigWig


def evaluate_track_predictions(
    single_task_data: dict,
    bin_size: int,
    ):
    
    """
    Calculates Pearson r and Jensen-Shannon divergence and extracts metadata for single prediction task,
    using a BIGWIG file and a single task data dictionary.
    
    Args:
        single_task_data (dict): The dictionary with predictions and metadata for a single task from the `prediction_tasks` 
                                 list of a predictions JSON.
        bin_size (int) : Returned from predictor
    Returns:
        correlation_details (dict): Dictionary containing 'pearson_r' (float, None, 0.0), Jensen-Shannon distance, and task metadata ("task_name", "task_type", "cell_type_actual")
    """
    # Extract metadata from single task data
    print(f"\n--- Extracting prediction_task metadata ---")
    task_name = single_task_data.get("name")
    task_type_actual = single_task_data.get("type_actual")
    cell_type_actual = single_task_data.get("cell_type_actual")
    predictions_dict =  single_task_data.get("predictions")

    
    scale_prediction_actual = single_task_data.get("scale_prediction_actual", None)
    trim_upstream_dict = single_task_data.get("trim_upstream", None)
    pearson_r_value = None # If there's an error, default to None
    
    if "error" in predictions_dict:
        print("No predictions were returned for this task -> Skipping evaluation calculation")
        return None

    # --- Data Validation and processing ---
    print("--- Validating data ---")
    if not isinstance(predictions_dict, dict):
        print(f"WARNING: 'predictions' in task: '{task_name}'\
            \nCell type: {cell_type_actual}\
            \nType: {task_type_actual}\
            \nis not a valid dictionary or is missing.")
    elif not predictions_dict:
        print(f"WARNING: 'predictions' dictionary is empty in task: '{task_name}'\
            \nCell type: {cell_type_actual}\
            \nType: {task_type_actual}")
    else:
        # Proceed with calculation if checks pass
        # Create DataFrame from Predictions
        print("--- Creating predictions_df ---")
        predictions_df = pd.DataFrame(list(predictions_dict.items()), columns=['region', 'Predicted_Value'])
        # check for NA in the predictions: a whole-cell None/NaN OR a NaN element inside a region's array
        na_mask = predictions_df['Predicted_Value'].apply(
            lambda v: v is None or np.isnan(np.asarray(v, dtype=float)).any()
        )
        na_rows = predictions_df[na_mask]
        if not na_rows.empty:
            print(f"NA values were found in the predictions ({len(na_rows)} region(s)), skipping evaluation")
            return None

        if trim_upstream_dict:
            trim_upstream_df = pd.DataFrame(list(trim_upstream_dict.items()), columns=['region', 'Trim Upstream Value'])
            predictions_df = pd.merge(predictions_df, trim_upstream_df, on = 'region', how = "left")
            predictions_df['Trim Upstream Value'] = predictions_df['Trim Upstream Value'].fillna(0)
            predictions_df['Trim Upstream Value'] = predictions_df['Trim Upstream Value'].astype(int)

            all_positive = (predictions_df['Trim Upstream Value'] >= 0).all()
            if all_positive == False:
                raise ValueError("Trim Upstream values are negative")
        else:
            predictions_df['Trim Upstream Value'] = 0
        

        if scale_prediction_actual == 'linear':
            print("Measured scale matches Predictor scale.")
        elif scale_prediction_actual is None or scale_prediction_actual == 'log':
            return None


        bigwig_file = pyBigWig.open(MEASURED_DATA_PATH)
        sequences = pd.read_parquet(EVALUATOR_INPUT_PATH)
        all_measurements = []
        for f in range(0, len(sequences)):
            chrom = sequences['chrom'].iloc[f]
            start = sequences['start'].iloc[f]
            end = sequences['end'].iloc[f]
            measurements = bigwig_file.values(chrom, start, end)
            #Mask missing values in the BigWig with 0
            measurements = np.nan_to_num(measurements, nan=0.0)
            all_measurements.append(measurements)

        sequences['measurements'] = all_measurements
        # Map region -> measurement so we align by name, not by position.
        # KeyError here (a prediction region absent from sequences) is intentional: fail loud.
        measurements_by_region = dict(zip(sequences['region'], sequences['measurements']))
        pearson_r_list = []
        jsd_list = []
        for i in range(0, len(predictions_df)):
            try:
                region_i = predictions_df['region'].iloc[i]
                measured = measurements_by_region[region_i]
                binned_predictions_current = predictions_df['Predicted_Value'].iloc[i]
                num_repeats = bin_size

                expanded_prediction_bp_level = [element for element in binned_predictions_current for _ in range(num_repeats)]

                #Once you expand the prediction to bp level, crop to match sequence length using trim upstream
                if predictions_df['Trim Upstream Value'].iloc[i] == 0:

                    length_expanded_predictions = len(expanded_prediction_bp_level)
                    length_measuremnts = len(measured)
                    trim_downstream = length_expanded_predictions - length_measuremnts

                    if trim_downstream == 0:
                        expanded_prediction_bp_level_trimmed = expanded_prediction_bp_level
                    else:
                        expanded_prediction_bp_level_trimmed = expanded_prediction_bp_level[:-trim_downstream]

                else:
                    trim_upstream = predictions_df['Trim Upstream Value'].iloc[i]
                    expanded_prediction_bp_level_trimmed_upstream = expanded_prediction_bp_level[trim_upstream:]
                    trim_downstream = len(expanded_prediction_bp_level_trimmed_upstream) - len(measured)

                    if trim_downstream == 0:
                        expanded_prediction_bp_level_trimmed = expanded_prediction_bp_level_trimmed_upstream
                    else:
                        expanded_prediction_bp_level_trimmed = expanded_prediction_bp_level_trimmed_upstream[:-trim_downstream]

                if len(expanded_prediction_bp_level_trimmed) == len(measured):
                    expanded_prediction_bp_level_trimmed = np.asarray(expanded_prediction_bp_level_trimmed, dtype=float)
                    measured = np.asarray(measured, dtype=float)

                    if expanded_prediction_bp_level_trimmed.std() == 0:
                        pearson_r_list.append(0.0)

                    else:
                        r, _ = pearsonr(expanded_prediction_bp_level_trimmed, measured)
                        pearson_r_list.append(0.0 if np.isnan(r) else float(r))

                    jsd = distance.jensenshannon(expanded_prediction_bp_level_trimmed, measured)
                    jsd_list.append(float(jsd))

                else:
                    print(f"WARNING: Length of predictions is not as expected so evaluation could not be completed")

            except ValueError as e:
                print(f"ValueError during Pearson correlation calculation for task: '{task_name}': {e}")

    if len(pearson_r_list) != 0 or len(jsd_list) != 0:
        correlation_details = {
            'task_name': task_name, 
            'task_type': task_type_actual,
            'cell_type_actual': cell_type_actual,
            'pearson_r': np.mean(pearson_r_list),
            'Jensen-Shannon distance': np.mean(jsd_list) if not np.any(np.isnan(jsd_list)) else "NaN"
        }
    else:
        print(f"WARNING: No evaluation metrics were calculated")
        correlation_details = {
            'task_name': task_name, 
            'task_type': task_type_actual,
            'cell_type_actual': cell_type_actual,
            'pearson_r': "NaN",
            'Jensen-Shannon distance': "NaN"
        }

    return correlation_details
    
def calculate_and_save_metrics(saved_predictions_path, output_dir, number_of_sequences):
    print(number_of_sequences)
    try:
        if os.path.exists(saved_predictions_path):
            print("----- Starting Evaluation Calculation and Saving as CSV -----")
            print(f"Using measured data from: {MEASURED_DATA_PATH}")
            print(f"Using predictions from: {saved_predictions_path}")
            print(f"Correlation metadata will be saved in {output_dir}")
            

            correlation_summary_filename = f"evaluation_summary_{EVALUATOR_NAME}.csv"
            correlation_summary_filepath = os.path.join(output_dir, correlation_summary_filename)
            
            # Initialize an empty list to get summary for all tasks
            all_task_correlation_results = []
            
            try:
              
                # Now load predictions
                with open(saved_predictions_path, 'r') as f:
                    predictions_file_content = json.load(f)

                # Extract Predictor Name
                predictor_name_base = predictions_file_content.get("predictor_name", "UnknownPredictor") # Resort UnknownPredictor if predictor name is not available
                predictor_name = predictor_name_base.replace(" ", "_").replace("/", "_") # Sanitize spaces and slashes for filenames, add underscores

                #extract model bin size
                bin_size = predictions_file_content.get("bin_size", 1)
                if (
                    "prediction_tasks" not in predictions_file_content or
                    # Also flag cases in case prediction_tasks key is returned empty
                    not predictions_file_content["prediction_tasks"] or
                    # And flag if any 'predictions' keys are empty
                    any(not key.get("predictions") for key in predictions_file_content["prediction_tasks"])
                ):
                    print("WARNING: 'prediction_tasks' key missing, empty, or one of the tasks has empty predictions.")
                else:
                    # Loop through each prediction_task from Predictor
                    # Calculate the correlation for each task separately
                    for task_index, single_task_data_dict in enumerate(predictions_file_content["prediction_tasks"]):
                        if not isinstance(single_task_data_dict, dict):
                            print(f"WARNING: Task item at index {task_index} is not a dictionary. Skipping!")
                            continue
                        
                        # Extract metadata from this task
                        task_type_actual = single_task_data_dict.get("type_actual")
                        predicted_cell_type = single_task_data_dict.get("cell_type_actual")
                        # We also want to extract the cell_type_requested to map it to measured_value_columns_map
                        requested_cell_type = single_task_data_dict.get("cell_type_requested")
                        
                        # Find the corresponding measured data column from the map
                        measured_col_for_task = 'measurements'
                        
                        print(f"\nProcessing task {task_index+1}")
                        prediction_task_data_nopredictions = [{k: v for k, v in single_task_data_dict.items() if k not in ("predictions", "trim_upstream")}]
                        # Call the correlation calculation function
                        task_correlation_dict = evaluate_track_predictions(
                            single_task_data=single_task_data_dict,
                            bin_size = bin_size # chromosome_column and chromosomes_to_filter_list can be added as arguments
                        )
                        
                        if task_correlation_dict:
                            pearson_r_value = task_correlation_dict.get('pearson_r')
                            jsd_value = task_correlation_dict.get('Jensen-Shannon distance')
                            # Get UTC timestamp for predictor_name
                            timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
                            # And append it to the predictor_name
                            predictor_identifier = f"{predictor_name_base}_{task_index}_{timestamp}" if predictor_name_base else f"UnknownPredictor_{task_index}_{timestamp}"
                            description = "DNase Track Request (K562)"
                            all_task_correlation_results.append({
                                "evaluator_name": EVALUATOR_NAME,
                                "description": description,
                                "predictor_name": predictor_name,
                                "time_stamp": timestamp,
                                'metric': 'pearson_r',
                                'value': str(pearson_r_value),
                                'prediction_task(s)_data': prediction_task_data_nopredictions,
                            })
                            all_task_correlation_results.append({
                                "evaluator_name": EVALUATOR_NAME,
                                "description": description,
                                "predictor_name": predictor_name,
                                "time_stamp": timestamp,
                                'metric': 'Jensen-Shannon distance',
                                'value': str(jsd_value),
                                'prediction_task(s)_data': prediction_task_data_nopredictions,
                            })

            except Exception as e:
                print(f"An error occurred during correlation calculation: {e}")
                
        # Once all the data is received, save them all into a summary CSV
            # print(all_task_correlation_results)
            if all_task_correlation_results:
                summary_df = pd.DataFrame(all_task_correlation_results)
                csv_file_exists: bool = os.path.isfile(correlation_summary_filepath)
                try:
                    summary_df.to_csv(correlation_summary_filepath, mode='a',
                                    sep='\t', header=(not csv_file_exists), index=False)
                    if csv_file_exists:
                        print("Appended to existing summary CSV file")
                    else:
                        print("Created a new summary CSV file")
                    print(f"Saved correlation summary to {correlation_summary_filepath}!")
                except IOError as e:
                    print("\nNo correlation resuls were saved!")

        else:
            print("Evaluator run did not complete successfully.")
            print(f"Predictions file not found in '{saved_predictions_path}'.")
            print("Skipping correlation calculation!")

    except Exception as e:
        print(f"An unexpected error occurred during evaluation calculations: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
