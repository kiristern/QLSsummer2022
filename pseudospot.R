# from https://github.com/ZJUFanLab/SpaTalk/blob/main/R/methods.R

#' @title Generate pseudo spot st_data
#'
#' @description Generate pseudo spot st_data with single-cell st_data
#' @param st_data A data.frame or matrix or dgCMatrix containing counts of spatial transcriptomics, each column representing a cell, each row representing a gene.
#' @param st_meta A data.frame containing coordinate of spatial transcriptomics with three columns, \code{'cell'}, \code{'x'}, \code{'y'}, and \code{celltype}.
#' @param x_min Min value of x axis.
#' @param x_res Resolution of x coordinate.
#' @param x_max Max value of x axis.
#' @param y_min Min value of y axis.
#' @param y_res Resolution of y coordinate.
#' @param y_max Max value of y axis.
#' @return A list of spot st_data and st_meta
#' @export
#' @import Matrix
#' @importFrom reshape2 dcast
#' @importFrom methods as

generate_spot <- function(st_data, st_meta, x_min, x_res, x_max, y_min, y_res, y_max) {
    if (is(st_data, "data.frame")) {
        st_data <- methods::as(as.matrix(st_data), "dgCMatrix")
    }
    if (is(st_data, "matrix")) {
        st_data <- methods::as(st_data, "dgCMatrix")
    }
    if (!is(st_data, "dgCMatrix")) {
        stop("st_data must be a data.frame or matrix or dgCMatrix!")
    }
    if (is(st_meta, "data.frame")) {
        if (!all(c("cell","x","y") %in% colnames(st_meta))) {
            stop("Please provide a correct st_meta data.frame! See demo_st_sc_meta()!")
        }
    } else {
        stop("st_meta must be a data.frame!")
    }
    x_range <- seq(from = x_min, to = x_max, by = x_res)
    y_range <- seq(from = y_min, to = y_max, by = y_res)
    # x and y resolution to construct spot data
    test_spot_plot <- data.frame(spot = "x", x = 0, y = 0, stringsAsFactors = F)
    for (i in 1:(length(x_range) - 1)) {
        x1 <- paste0("x", rep(i, length(y_range) - 1))
        test_spot_plot1 <- data.frame(spot = x1, x = x_range[i] + x_res, y = 0, stringsAsFactors = F)
        for (j in 1:(length(y_range) - 1)) {
            y1 <- paste0(test_spot_plot1$spot[j], "_", "y", j)
            test_spot_plot1$spot[j] <- y1
            test_spot_plot1$y[j] <- y_range[j] + y_res
        }
        test_spot_plot <- rbind(test_spot_plot, test_spot_plot1)
    }
    test_spot_plot <- test_spot_plot[-1, ]
    for (i in 1:nrow(st_meta)) {
        test_data_x <- max(which(x_range <= st_meta$x[i]))
        test_data_y <- max(which(y_range <= st_meta$y[i]))
        st_meta$spot[i] <- paste0("x", test_data_x, "_", "y", test_data_y)
    }
    test_spot_meta <- as.data.frame(table(st_meta$spot), stringsAsFactors = F)
    test_spot_meta1 <- reshape2::dcast(data = st_meta[, c("spot", "celltype")], formula = spot ~ celltype)
    test_spot_meta <- cbind(test_spot_meta, test_spot_meta1[, -1])
    # test_spot_data -- sum
    test_spot_data <- list()
    for (i in 1:nrow(test_spot_meta)) {
        test_spot_cell <- st_data[, st_meta[st_meta$spot == test_spot_meta$Var1[i], ]$cell]
        if (is(test_spot_cell, "dgCMatrix")) {
            test_spot_sum <- rowSums(test_spot_cell)
        } else {
            test_spot_sum <- test_spot_cell
        }
        test_spot_data[[i]] <- test_spot_sum
        names(test_spot_data)[i] <- test_spot_meta$Var1[i]
    }
    test_spot_data <- as.data.frame(test_spot_data, stringsAsFactors = F)
    test_spot_data <- as(as.matrix(test_spot_data), "dgCMatrix")
    # generate x and y
    rownames(test_spot_plot) <- test_spot_plot$spot
    test_spot_plot <- test_spot_plot[test_spot_meta$Var1, ]
    test_spot_real <- test_spot_meta
    colnames(test_spot_real)[c(1, 2)] <- c("spot", "cell_real")
    test_spot_meta <- test_spot_plot
    test_spot_meta <- cbind(test_spot_meta, test_spot_real[, -1])
    return(list(st_data = test_spot_data, st_meta = test_spot_meta))
}

