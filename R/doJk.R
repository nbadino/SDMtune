#' Jackknife Test
#'
#' Run the Jackknife test for variable importance removing one variable at time.
#'
#' @param model \linkS4class{SDMmodel} or \linkS4class{SDMmodelCV} object.
#' @param metric character. The metric used to evaluate the models, possible
#' values are: "auc", "tss" and "aicc".
#' @param variables vector. Variables used for the test, if not provided it
#' takes all the variables used to train the model.
#' @param test \linkS4class{SWD}. If provided it reports the result also for the
#' testing dataset. Not used for **aicc** and \linkS4class{SDMmodelCV}.
#' @param with_only logical. If `TRUE` it runs the test also for each variable
#' in isolation.
#' @param env \link[terra]{rast} containing the environmental variables, used
#' only with "aicc".
#' @param return_models logical. If `TRUE` returns all the models together with
#' the test result.
#' @param progress logical If `TRUE` shows a progress bar.
#' @param parallel logical. If `TRUE` trains models in parallel across variables.
#' @param ncores integer. Number of cores to use for parallel training when
#' `parallel = TRUE`. Default is all available cores minus one.
#'
#' @return A data frame with the test results. If `return_model = TRUE` it
#' returns a list containing the test results together with the models.
#' @export
#'
#' @author Sergio Vignali
#'
#' @examples
#' \donttest{# Acquire environmental variables
#' files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
#'                     pattern = "grd",
#'                     full.names = TRUE)
#'
#' predictors <- terra::rast(files)
#'
#' # Prepare presence and background locations
#' p_coords <- virtualSp$presence
#' bg_coords <- virtualSp$background
#'
#' # Create SWD object
#' data <- prepareSWD(species = "Virtual species",
#'                    p = p_coords,
#'                    a = bg_coords,
#'                    env = predictors,
#'                    categorical = "biome")
#'
#' # Split presence locations in training (80%) and testing (20%) datasets
#' datasets <- trainValTest(data,
#'                          test = 0.2,
#'                          only_presence = TRUE)
#' train <- datasets[[1]]
#' test <- datasets[[2]]
#'
#' # Train a model
#' model <- train(method = "Maxnet",
#'                data = train,
#'                fc = "lq")
#'
#' # Execute the Jackknife test only for the environmental variables "bio1" and
#' # "bio12", using the metric AUC
#' doJk(model,
#'      metric = "auc",
#'     variables = c("bio1", "bio12"),
#'     test = test)
#'
#' # The same without testing dataset
#' doJk(model,
#'      metric = "auc",
#'      variables = c("bio1", "bio12"))
#'
#' # Execute the Jackknife test only for the environmental variables "bio1" and
#' # "bio12", using the metric TSS but without running the test for one single
#' # variable
#' doJk(model,
#'      metric = "tss",
#'      variables = c("bio1", "bio12"),
#'      test = test,
#'      with_only = FALSE)
#'
#' # Execute the Jackknife test only for the environmental variables "bio1" and
#' # "bio12", using the metric AICc but without running the test for one single
#' # variable
#' doJk(model,
#'      metric = "aicc",
#'      variables = c("bio1", "bio12"),
#'      with_only = FALSE,
#'      env = predictors)
#'
#' # Execute the Jackknife test for all the environmental variables using the
#' # metric AUC and returning all the trained models
#' jk <- doJk(model,
#'            metric = "auc",
#'            test = test,
#'            return_models = TRUE)
#'
#' jk$results
#' jk$models_without
#' jk$models_withonly}
doJk <- function(model,
                 metric,
                 variables = NULL,
                 test = NULL,
                 with_only = TRUE,
                 env = NULL,
                 return_models = FALSE,
                 progress = TRUE,
                 parallel = FALSE,
                 ncores = parallel::detectCores() - 1) {

  metric <- match.arg(metric, c("auc", "tss", "aicc"))

  .check_args(model, metric = metric, test = test, env = env, is_do_jk = TRUE)

  if (is.null(variables))
    variables <- colnames(model@data@data)

  n <- length(variables)

  if (with_only) {
    tot <- n * 2
  } else {
    tot <- n
  }

  if (!is.null(env)) {
    # TODO: Remove with version 2.0.0
    if (inherits(env, "Raster")) {
      .raster_error("rast")
    }

    if (!inherits(env, "SpatRaster"))
      cli::cli_abort(c(
        "!" = "{.var env} must be a {.cls SpatRaster} object",
        "x" = "You have supplied a {.cls {class(env)}} instead."
      ))
  }

  if (progress)
    cli::cli_progress_bar(
      name = "Jk Test",
      type = "iterator",
      format = "{cli::pb_name} {cli::pb_bar} {cli::pb_percent} | \\
                ETA: {cli::pb_eta} - {cli::pb_elapsed_clock}",
      total = tot,
      clear = FALSE
    )

  models_without <- vector("list", length = n)
  models_withonly <- vector("list", length = n)

  res <- matrix(nrow = n, ncol = 5)

  if (metric == "auc") {
    labels <- c("Variable", "Train_AUC_without", "Train_AUC_withonly",
                "Test_AUC_without", "Test_AUC_withonly")
  } else if (metric == "tss") {
    labels <- c("Variable", "Train_TSS_without", "Train_TSS_withonly",
                "Test_TSS_without", "Test_TSS_withonly")
  } else {
    labels <- c("Variable", "AICc_without", "AICc_withonly", "-", "-")
  }

  old_model <- model

  if (inherits(model, "SDMmodelCV"))
    t <- TRUE

  do_parallel <- parallel && n > 1

  if (do_parallel) {
    if (progress)
      cli::cli_alert_info("Running Jackknife test in parallel on {min(ncores, n)} core{?s}...")

    if (.Platform$OS.type == "windows") {
      ncores <- min(ncores, n)
      cl <- parallel::makePSOCKcluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("old_model", "variables", "metric", "test",
                                    "env", "with_only"),
                              envir = environment())
      parallel::clusterEvalQ(cl, library(SDMtune))
      results <- parallel::clusterApplyLB(cl, seq_len(n), function(i) {
        var <- variables[i]
        data_wo <- old_model@data
        data_wo@data[var] <- NULL
        t_wo <- NULL
        if (metric != "aicc" & !inherits(old_model, "SDMmodelCV") & !is.null(test)) {
          t_wo <- test
          t_wo@data[var] <- NULL
        }
        settings <- list("data" = data_wo)
        jk_model_wo <- SDMtune:::.create_model_from_settings(old_model, settings)
        train_wo <- SDMtune:::.get_metric(metric, jk_model_wo, env = env)
        test_wo <- NA_real_
        if (metric != "aicc" & !is.null(test))
          test_wo <- SDMtune:::.get_metric(metric, jk_model_wo, test = t_wo)
        train_wo_only <- NA_real_
        test_wo_only <- NA_real_
        jk_model_wo_only <- NULL
        if (with_only) {
          data_wo_only <- old_model@data
          data_wo_only@data <- data_wo_only@data[var]
          t_wo_only <- NULL
          if (metric != "aicc" & !inherits(old_model, "SDMmodelCV") & !is.null(test)) {
            t_wo_only <- test
            t_wo_only@data <- t_wo_only@data[var]
          }
          settings_only <- list("data" = data_wo_only)
          jk_model_wo_only <- SDMtune:::.create_model_from_settings(old_model, settings_only)
          train_wo_only <- SDMtune:::.get_metric(metric, jk_model_wo_only, env = env)
          if (metric != "aicc" & !is.null(test))
            test_wo_only <- SDMtune:::.get_metric(metric, jk_model_wo_only, test = t_wo_only)
        }
        list(train_wo = train_wo, train_wo_only = train_wo_only,
             test_wo = test_wo, test_wo_only = test_wo_only,
             model_wo = jk_model_wo, model_wo_only = jk_model_wo_only)
      })
    } else {
      ncores <- min(ncores, n)
      results <- parallel::mclapply(seq_len(n), function(i) {
        Sys.setenv(OMP_NUM_THREADS = "1")
        var <- variables[i]
        data_wo <- old_model@data
        data_wo@data[var] <- NULL
        t_wo <- NULL
        if (metric != "aicc" & !inherits(old_model, "SDMmodelCV") & !is.null(test)) {
          t_wo <- test
          t_wo@data[var] <- NULL
        }
        settings <- list("data" = data_wo)
        jk_model_wo <- .create_model_from_settings(old_model, settings)
        train_wo <- .get_metric(metric, jk_model_wo, env = env)
        test_wo <- NA_real_
        if (metric != "aicc" & !is.null(test))
          test_wo <- .get_metric(metric, jk_model_wo, test = t_wo)
        train_wo_only <- NA_real_
        test_wo_only <- NA_real_
        jk_model_wo_only <- NULL
        if (with_only) {
          data_wo_only <- old_model@data
          data_wo_only@data <- data_wo_only@data[var]
          t_wo_only <- NULL
          if (metric != "aicc" & !inherits(old_model, "SDMmodelCV") & !is.null(test)) {
            t_wo_only <- test
            t_wo_only@data <- t_wo_only@data[var]
          }
          settings_only <- list("data" = data_wo_only)
          jk_model_wo_only <- .create_model_from_settings(old_model, settings_only)
          train_wo_only <- .get_metric(metric, jk_model_wo_only, env = env)
          if (metric != "aicc" & !is.null(test))
            test_wo_only <- .get_metric(metric, jk_model_wo_only, test = t_wo_only)
        }
        list(train_wo = train_wo, train_wo_only = train_wo_only,
             test_wo = test_wo, test_wo_only = test_wo_only,
             model_wo = jk_model_wo, model_wo_only = jk_model_wo_only)
      }, mc.cores = ncores)
    }

    for (i in seq_len(n)) {
      r <- results[[i]]
      res[i, 2] <- r$train_wo
      res[i, 3] <- r$train_wo_only
      res[i, 4] <- r$test_wo
      res[i, 5] <- r$test_wo_only
      models_without[[i]] <- r$model_wo
      if (with_only)
        models_withonly[[i]] <- r$model_wo_only
      if (progress) {
        cli::cli_progress_update()
        if (with_only) cli::cli_progress_update()
      }
    }

  } else {

    for (i in 1:n) {
      data <- old_model@data
      data@data[variables[i]] <- NULL

      if (metric != "aicc" & !inherits(model, "SDMmodelCV") & !is.null(test)) {
        t <- test
        t@data[variables[i]] <- NULL
      }

      settings <- list("data" = data)

      jk_model <- .create_model_from_settings(model, settings)

      res[i, 2] <- .get_metric(metric, jk_model, env = env)

      if (metric != "aicc" & !is.null(test))
        res[i, 4] <- .get_metric(metric, jk_model, test = t)

      models_without[[i]] <- jk_model

      if (progress)
        cli::cli_progress_update()

      if (with_only) {
        data <- old_model@data
        data@data <- data@data[variables[i]]

        if (metric != "aicc" & !inherits(model, "SDMmodelCV") & !is.null(test)) {
          t <- test
          t@data <- t@data[variables[i]]
        }

        settings <- list("data" = data)

        jk_model <- .create_model_from_settings(model, settings)

        res[i, 3] <- .get_metric(metric, jk_model, env = env)

        if (metric != "aicc" & !is.null(test))
          res[i, 5] <- .get_metric(metric, jk_model, test = t)

        models_withonly[[i]] <- jk_model

        if (progress)
          cli::cli_progress_update()
      }
    }
  }

  jk_test <- as.data.frame(res, stringAsFactor = FALSE)
  colnames(jk_test) <- labels
  jk_test[1] <- variables

  jk_test <- Filter(function(x) !all(is.na(x)), jk_test)

  if (return_models) {

    if (with_only) {
      output <- list(results = jk_test, models_without =  models_without,
                     models_withonly = models_withonly)
    } else {
      output <- list(results = jk_test, models_without = models_without)
    }

  } else {
    output <- jk_test
  }

  output
}
