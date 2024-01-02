#' @export
rhoControlUI <- function(p, input) {
    sp <- p@species_params[input$sp, ]
    tagList(
        tags$h3(tags$a(id = "rho"), "rho"),
        # sliderInput("rho_detritus", "rho_detritus", value = sp$rho_detritus,
        #             min = 0,
        #             max = signif(ifelse(sp$rho_detritus > 0,
        #                                 sp$rho_detritus * 2,
        #                                 0.001), 2)),
        # sliderInput("n_detritus", "n_detritus", value = sp$n_detritus,
        #             min = -.5,
        #             max = .75, step = .05),
        sliderInput("rho_algae", "rho_algae", value = sp$rho_algae,
                    min = 0,
                    max = signif(ifelse(sp$rho_algae > 0,
                                        sp$rho_algae * 2,
                                        0.001), 2))
    )
}

#' @export
rhoControl <- function(input, output, session, params, params_old, flags, ...) {
    observeEvent(
        list(input$rho_detritus, input$rho_algae),
        {
            p <- params()
            sp <- input$sp
            if (!identical(sp, flags$sp_old_rho)) {
                flags$sp_old_rho <- sp
                return()
            }
            # Update slider min/max so that they are a fixed proportion of the
            # parameter value
            # updateSliderInput(session, "rho_detritus",
            #                   min = 0,
            #                   max = signif(ifelse(input$rho_detritus > 0,
            #                                       input$rho_detritus * 2,
            #                                       0.001), 2))
            updateSliderInput(session, "rho_algae",
                              min = 0,
                              max = signif(ifelse(input$rho_algae > 0,
                                                  input$rho_algae * 2,
                                                  0.001), 2))
            
            # p@species_params[sp, "rho_detritus"]   <- input$rho_detritus
            # p@species_params[sp, "n_detritus"]   <- input$n_detritus
            p@species_params[sp, "rho_algae"]   <- input$rho_algae
            p <- setRho(p)
            mizerExperimental:::tuneParams_update_species(sp, p, params, params_old)
        },
        ignoreInit = TRUE)
}