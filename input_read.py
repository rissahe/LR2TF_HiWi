
#' Loading old dorothea package regulons (default if no regulon is provided)
#'
#' @param arguments_list list of user defined arguments
#' @return val_arguments validated list of arguments
#' @export

def validate_input_arguments (arguments_list,**kwargs):
    if out_path in arguments_list is None:
        print("Please provide an output path")
    else: 
        if 