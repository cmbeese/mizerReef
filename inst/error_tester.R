param_tester <- function(params, 
                         species, 
                         parameter, 
                         param_list) {
    
    library(tryCatchLog)
    
    # Initialize a dataframe to store errors
    errors_df <- data.frame(Parameter = character(), 
                            Error = character(), 
                            stringsAsFactors = FALSE)
    
    for (i in 1:length(param_list)) {
        tryCatch({
            # Attempt to run the function with the current parameter
            species_params(params)[species, parameter] <- param_list[i]
            matchReefGrowth(params)
        }, error = function(e) {
            # If an error occurs, store the parameter and error message
            errors_df <- rbind(errors_df, 
                               data.frame(Parameter = as.character(parameter), 
                                          ErrorMessage = e$message))
        })
    }
    return(errors_df)
}
