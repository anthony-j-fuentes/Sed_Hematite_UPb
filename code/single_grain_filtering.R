
library(dplyr)
library(stringr)
library(IsoplotR)
library(ggplot2)

setwd("~/GitHub/Sed_Hematite_UPb")


##### Functions ####
splt_data <- function(df, sample, prop) {
  # Create a new column 'sample_grp' by stripping the trailing underscore and digits
  df <- df %>% mutate(sample_grp = str_remove(X, "_\\d+$"))
  
  if (prop) {
    samp <- df %>%
      filter(sample_grp == sample) %>%
      select(
        `X238U.206Pb` = "Final.U238.Pb206_mean",
        `X238U.206Pb.2SE.prop.` = "Final.U238.Pb206_2SE.prop.",
        `X207Pb.206Pb` = "Final.Pb207.Pb206_mean",
        `X207Pb.206Pb.2SE.prop.` = "Final.Pb207.Pb206_2SE.prop.",
        `Err.Corr.7.6.v.8.6` = "rho.207Pb.206Pb.v.238U.206Pb"
      ) %>%
      #distinct() %>%
      na.omit()
    return(read.data(samp, method = "U-Pb", format = 2, ierr = 2))
  } else {
    samp <- df %>%
      filter(sample_grp == sample) %>%
      select(
        `X238U.206Pb` = "Final.U238.Pb206_mean",
        `X238U.206Pb.2SE.int.` = "Final.U238.Pb206_2SE.int.",
        `X207Pb.206Pb` = "Final.Pb207.Pb206_mean",
        `X207Pb.206Pb.2SE.int.` = "Final.Pb207.Pb206_2SE.int.",
        `Err.Corr.7.6.v.8.6` = "rho.207Pb.206Pb.v.238U.206Pb"
      ) %>%
      #distinct() %>%
      na.omit()
    return(read.data(samp, method = "U-Pb", format = 2, ierr = 2))
  }
}

## This function produces summary data prior to filtering 
get_spread_report <- function(df, name_filter) {
  df <- df %>%
    mutate(sample_grp = str_remove(X, "_\\d+$")) %>%
    filter(str_detect(X, name_filter))
  
  report_list <- list()
  
  for (sg in unique(df$sample_grp)) {
    df_sub <- df %>% filter(sample_grp == sg)
    if (nrow(df_sub) < 2) next
    
    # Spread Calculations (required for audit_thrs)
    p_low_u  <- df_sub[which.min(df_sub$Final.U238.Pb206_mean), ]
    p_high_u <- df_sub[which.max(df_sub$Final.U238.Pb206_mean), ]
    mean_x   <- mean(df_sub$Final.U238.Pb206_mean, na.rm = TRUE)
    mean_y   <- mean(df_sub$Final.Pb207.Pb206_mean, na.rm = TRUE)
    
    abs_u_spread  <- abs(p_high_u$Final.U238.Pb206_mean - p_low_u$Final.U238.Pb206_mean)
    abs_pb_spread <- abs(p_high_u$Final.Pb207.Pb206_mean - p_low_u$Final.Pb207.Pb206_mean)
    
    age_val <- mswd_val <- ci_val <- pb_int <- NaN
    
    tryCatch({
      grain_data <- splt_data(df, sg, prop = FALSE)
      # oerr=2 ensures the CI is returned directly
      indiv_age  <- age(grain_data, method = "U238-Pb206", 
                        type = 3, common.Pb = 0, isochron = TRUE, oerr = 2)
      
      if (!is.null(indiv_age$par)) {
        age_val  <- as.numeric(indiv_age$par['t'])
        ci_val   <- as.numeric(indiv_age$err[2, 't'])
        pb_int   <- as.numeric(indiv_age$par['76i'])
        mswd_val <- as.numeric(indiv_age$mswd)
      }
    }, error = function(e) {})
    
    report_list[[sg]] <- data.frame(
      Sample_Name = sg,
      U_Pb_Abs_Spread = round(abs_u_spread, 4),
      Pb_Pb_Abs_Spread = round(abs_pb_spread, 4),
      U_Pb_Pct_Spread  = round((abs_u_spread/mean_x)*100, 2),
      Pb_Pb_Pct_Spread = round((abs_pb_spread/mean_y)*100, 2),
      Valid_Slope = (p_high_u$Final.Pb207.Pb206_mean < p_low_u$Final.Pb207.Pb206_mean),
      n_points = nrow(df_sub),
      age = age_val, mswd = mswd_val, ci = ci_val, pb_intercept = pb_int,
      stringsAsFactors = FALSE
    )
  }
  return(bind_rows(report_list))
}


plot_spread_v_error <- function(report_df, name, type = "Pb", scale = "relative") {
  library(ggplot2)
  
  # 1. Prepare data: Calculate Relative Error
  plot_data <- report_df %>%
    filter(!is.na(age), !is.na(ci), age > 0) %>%
    mutate(rel_error_pct = (ci / age))
  

  if (type == "U") {
    prefix <- "U/Pb"
    if (scale == "relative") {
      x_col <- "U_Pb_Pct_Spread"
      x_label <- "Normalized U/Pb Spread (%)"
    } else {
      x_col <- "U_Pb_Abs_Spread"
      x_label <- "Absolute U/Pb Spread "
    }
  } else if (type == "Pb") {
    prefix <- "Pb/Pb"
    if (scale == "relative") {
      x_col <- "Pb_Pb_Pct_Spread"
      x_label <- "Normalized Pb/Pb Spread (%)"
    } else {
      x_col <- "Pb_Pb_Abs_Spread"
      x_label <- "Absolute Pb/Pb Spread "
    }
  } else {
    stop("Type must be 'U' or 'Pb'")
  }
  

  p <- ggplot(plot_data, aes(x = .data[[x_col]], y = rel_error_pct)) +
    geom_point(aes(size = n_points, color = Valid_Slope), alpha = 0.7) +
    scale_color_manual(values = c("TRUE" = "#2c7bb6", "FALSE" = "#d7191c")) +
    theme_bw() +
    labs(
      title = paste(name, "Isochron Precision vs.", prefix, scale, "Spread"),
      subtitle = paste("Comparison of", x_label, "to Relative Age Error"),
      x = x_label,
      y = "Relative Age Error (95% CI / Age)",
      size = "n points",
      color = "Slope"
    )
  
  return(p)
}

audit_thrs <- function(report_df, start, end, step, 
                       type = "Pb", scale = "relative", 
                       precision_limit = 0.5) {
  library(dplyr)
  
  if (type == "U") {
    x_col <- if(scale == "relative") "U_Pb_Pct_Spread" else "U_Pb_Abs_Spread"
  } else {
    x_col <- if(scale == "relative") "Pb_Pb_Pct_Spread" else "Pb_Pb_Abs_Spread"
  }
  
  threshold_steps <- seq(from = start, to = end, by = step)
  

  results <- lapply(threshold_steps, function(t) {
    

    passed_data <- report_df %>%
      filter(.data[[x_col]] >= t) %>%
      filter(!is.na(age), !is.na(ci))
    
    total_grains <- nrow(passed_data)

    if (total_grains == 0) {
      return(data.frame(
        threshold = t,
        total_grains = 0,
        precise_grains = 0,
        success_proportion = NaN,
        stringsAsFactors = FALSE
      ))
    }
    

    precise_count <- sum((passed_data$ci / passed_data$age) <= precision_limit, na.rm = TRUE)
    
    data.frame(
      threshold = t,
      total_grains = total_grains,
      precise_grains = precise_count,
      success_proportion = precise_count / total_grains,
      stringsAsFactors = FALSE
    )
  })
  
  return(bind_rows(results))
}



plot_audit_results <- function(audit_df, type = "Pb", scale = "relative", 
                               hex_color = "#008080") {
  library(ggplot2)
  

  if (type == "Pb") {
    iso_math <- bquote(""^207*Pb/""^206*Pb)
  } else {
    iso_math <- bquote(""^238*U/""^206*Pb)
  }
  
  scale_txt <- if(scale == "relative") "Relative (%)" else "Absolute"
  x_lab_final <- bquote("Threshold (" * .(scale_txt) ~ .(iso_math) * ")")

  p <- ggplot(audit_df, aes(x = threshold, y = success_proportion)) +

    geom_line(aes(group = 1), color = hex_color, linewidth = 1) +
    

    geom_point(shape = 23, 
               size = 4, 
               fill = hex_color, 
               color = "black", 
               stroke = 0.8) +
    
    scale_y_continuous(
      limits = c(0, 1.05), 
      breaks = seq(0, 1, 0.1),
      expand = c(0, 0.02)
    ) +
    
    theme_bw() +
    labs(
      title = "Isochron Success Rate Analysis",
      subtitle = "Proportion of samples meeting precision goal (CI/Age < 0.5)",
      x = x_lab_final,
      y = "Success Proportion"
    ) +
    theme(
      axis.title = element_text(size = 11),
      plot.title = element_text(face = "bold")
    )
  
  return(p)
}


sample_spread_filter <- function(df, thrs,name) {
  df <- df %>% mutate(sample_grp = str_remove(X, "_\\d+$"))
  df<-df%>% filter(str_detect(X,name))
  sample_groups <- unique(df$sample_grp)
  filt_df <- data.frame()
  
  for (sg in sample_groups) {
    df_sub <- df %>% filter(sample_grp == sg)
    if (nrow(df_sub) <= 2) {
      next}
    
    spread <- (max(df_sub$Final.U238.Pb206_mean, na.rm = TRUE) - min(df_sub$Final.U238.Pb206_mean, na.rm = TRUE))
    
    if (spread >= thrs) {
      filt_df <- rbind(filt_df, df_sub)
    }
  }
  
  return(filt_df)
}



filter_by_spread <- function(df, thrs, name_filter, type = "U", scale = "relative") {
  
  df <- df %>%
    mutate(sample_grp = str_remove(X, "_\\d+$")) %>%
    filter(str_detect(X, name_filter))
  
  filt_list <- list()
  
  for (sg in unique(df$sample_grp)) {
    df_sub <- df %>% filter(sample_grp == sg)
    

    if (nrow(df_sub) <= 2) next

    p_low_u  <- df_sub[which.min(df_sub$Final.U238.Pb206_mean), ]
    p_high_u <- df_sub[which.max(df_sub$Final.U238.Pb206_mean), ]
    

    if (type == "U") {
      vals <- df_sub$Final.U238.Pb206_mean
    } else if (type == "Pb") {
      vals <- df_sub$Final.Pb207.Pb206_mean
    } else {
      stop("Type must be 'U' or 'Pb'")
    }
    

    abs_diff <- max(vals) - min(vals)
    
    if (scale == "relative") {
      actual_spread <- (abs_diff / mean(vals)) * 100
    } else if (scale == "absolute") {
      actual_spread <- abs_diff
    } else {
      stop("Scale must be 'relative' or 'absolute'")
    }
    

    if (p_high_u$Final.Pb207.Pb206_mean < p_low_u$Final.Pb207.Pb206_mean) {
      if (actual_spread >= thrs) {
        filt_list[[sg]] <- df_sub
      }
    }
  }
  
  return(bind_rows(filt_list))
}


mswd_filt <- function(df, prop, thrs, age_prp) {
  sample_groups <- unique(df$sample_grp)
  filt_mswd <- data.frame()
  
  for (sg in sample_groups) {
    indiv_age <- tryCatch({
      grain_data <- splt_data(df, sg, prop)

      age(grain_data, method = "U238-Pb206", type = 3, 
          common.Pb = 0, isochron = TRUE, oerr = 2)
    }, error = function(e) { return(NULL) })
    
    if (is.null(indiv_age) || is.null(indiv_age$par)) next
    

    age_val  <- as.numeric(indiv_age$par['t'])

    ci_val   <- as.numeric(indiv_age$err[2, 't']) 
    mswd_val <- as.numeric(indiv_age$mswd[1])
    

    if (is.na(mswd_val) || mswd_val > thrs) next
    if (is.na(age_val) || age_val > 4000) next
    if (is.na(ci_val) || ci_val > (age_prp * age_val)) next
    
    filt_mswd <- rbind(filt_mswd, data.frame(
      sample = sg, age = age_val, mswd = mswd_val,
      ci = ci_val, stringsAsFactors = FALSE
    ))
  }
  return(filt_mswd)
}

plot_wetherill_limited<-function(data_obj, label, xlim, ylim, save_pdf, filename) {
  concordia(data_obj, type = 1,anchor = c(2,0),
            ellipse.fill = alpha('#de2d26', 0.6), ticks = 10,
            show.numbers = FALSE, concordia.col = 'black',
            show.age = 2, oerr = 2, common.Pb = 0,
            cex.lab = 2, cex.axis = 1.5,
            xlim = xlim, ylim = ylim)
  usr <- par("usr")
  text(usr[2], usr[4], label, adj = c(1, 1), cex = 1.2)
  if (save_pdf)
    dev.copy2pdf(file = filename, width = 8, height = 6)
}



plot_by_sample <- function(df, prop, save_pdf = TRUE) {
  df <- df %>% mutate(sample_grp = str_remove(X, "_\\d+$"))
  sample_groups <- unique(df$sample_grp)
  cat("Plotting", length(sample_groups), "sample groups:\n")
  print(sample_groups)
  
  for (sg in sample_groups) {
    grain_data <- splt_data(df, sg, prop)
    
    tryCatch({
      plot_concordia_default(grain_data, paste("Grain:", sg),
                             save_pdf,
                             file.path(out_dir, paste0("concordia_", sg, ".pdf")))
    }, error = function(e) {
      cat("Error with sample", sg, ":", e$message, "\n")
    })
    
    tryCatch({
      # Using your preferred xlim/ylim values for the limited plot
      plot_concordia_limited(grain_data, paste("Grain:", sg),
                             xlim = c(0.1, 9), ylim = c(0.05, 0.9),
                             save_pdf,
                             file.path(out_dir, paste0("concordia_", sg, "_zoomout.pdf")))
    }, error = function(e) {
      cat("Error with sample", sg, ":", e$message, "\n")
    })
  }
}






plot_wetherill_sample <- function(df, prop, save_pdf = TRUE) {
  df <- df %>% mutate(sample_grp = str_remove(X, "_\\d+$"))
  sample_groups <- unique(df$sample_grp)
  cat("Plotting", length(sample_groups), "sample groups:\n")
  print(sample_groups)
  
  for (sg in sample_groups) {
    grain_data <- splt_data(df, sg, prop)
    
    tryCatch({
      plot_wetherill_default(grain_data, paste("Grain:", sg),
                             save_pdf,
                             file.path(out_dir, paste0("concordia_", sg, ".pdf")))
    }, error = function(e) {
      cat("Error with sample", sg, ":", e$message, "\n")
    })
    
    tryCatch({
      # Using your preferred xlim/ylim values for the limited plot
      plot_wetherill_limited(grain_data, paste("Grain:", sg),
                             xlim = c(0.1, 9), ylim = c(0.05, 0.9),
                             save_pdf,
                             file.path(out_dir, paste0("concordia_", sg, "_zoomout.pdf")))
    }, error = function(e) {
      cat("Error with sample", sg, ":", e$message, "\n")
    })
  }
}

plot_grains_wetherill<- function(df, grain,prop, hide,pts,zoom){
  
  
  
  if(hide==FALSE){
    grain_data<-splt_data(df,grain,prop)
    
    concordia(grain_data,type=1,anchor=c(2,0),ellipse.fill = alpha('#de2d26', 0.6), ticks = 10,
              show.numbers = TRUE,show.age = 2, oerr = 2, common.Pb = 0,cex.lab=(2),cex.axis=(1.5))
    usr <- par("usr")
    text(usr[2], usr[4], grain, adj = c(1, 1), cex = 1.2)
    
  }
  
  if(hide==TRUE){
    grain_data<-splt_data(df,grain,prop)
    print(rownames(grain_data))
    concordia(grain_data,type=1,anchor=c(2,0),ellipse.fill = alpha('#de2d26', 0.6), ticks = 10,
              show.numbers = TRUE,show.age = 2, oerr = 2, common.Pb = 0,omit=pts,cex.lab=(2),cex.axis=(1.5))
    
    usr <- par("usr")
    text(usr[2], usr[4], grain, adj = c(1, 1), cex = 1.2)
    
    
  }
  
  
}


# This function is for plotting each grain individually that I have screened using the spread filter above 

plot_grains<- function(df, grain,prop, hide,pts,zoom){
  

  
  if(hide==FALSE){
    grain_data<-splt_data(df,grain,prop)
   
    concordia(grain_data,type=2,ellipse.fill = alpha('#de2d26', 0.6), ticks = 10,
              show.numbers = TRUE,show.age = 2, oerr = 3, common.Pb = 0,cex.lab=(2),cex.axis=(1.5))
    usr <- par("usr")
    text(usr[2], usr[4], grain, adj = c(1, 1), cex = 1.2)
    
  }
    
  if(hide==TRUE){
    grain_data<-splt_data(df,grain,prop)
    print(rownames(grain_data))
    concordia(grain_data,type=2,ellipse.fill = alpha('#de2d26', 0.6), ticks = 10,
              show.numbers = TRUE,show.age = 2, oerr = 3, common.Pb = 0,omit=pts,
              cex.lab=(2),cex.axis=(1.5))
    
    usr <- par("usr")
    text(usr[2], usr[4], grain, adj = c(1, 1), cex = 1.2)
    
  
  }
  

}




remove_pts <- function(df, pts) {
  filtered_df <- df  
  
  for (i in 1:nrow(pts)) {
    sample <- pts$grain[i]
    omit_spts <- as.character(pts$spots[i])  # Fix: grab only the current row's spot list
    
    rmv_spt <- as.numeric(unlist(strsplit(omit_spts, ",")))
    spts_to_remove <- paste0(sample, "_", rmv_spt)
    print(spts_to_remove)
    # Remove matching sample spots
    filtered_df <- filtered_df[!filtered_df$X %in% spts_to_remove, ]
  }
  
  return(filtered_df)
}

### This function is for after screening for points to exclude and allows you to make zoomed plots for those grains 
## and specify points to remove

make_grain_plots<- function(df, grain,prop,pts,save){
  

    grain_data<-splt_data(df,grain,prop)
    concordia(grain_data,type=2,ellipse.fill = alpha('#de2d26', 0.6), ticks = 10,concordia.col = 'black',
              show.numbers = FALSE,show.age = 2, oerr = 2, common.Pb = 0,omit=pts,cex.lab=(2),cex.axis=(1.5),
              xlim = c(0.1, 9), ylim = c(0.05, 0.9))
    
    usr <- par("usr")
    text(usr[2], usr[4], grain, adj = c(1.2, 1.2), cex = 1.2)
    if (save==TRUE){
      filename=file.path(out_dir, paste0( grain,"full", ".pdf"))
      dev.copy2pdf(file = filename, width = 8, height = 6)
    }
    concordia(grain_data,type=2,ellipse.fill = alpha('#de2d26', 0.6), ticks = 10,concordia.col = 'black',
              show.numbers = FALSE,show.age = 2, oerr = 2, common.Pb = 0,omit=pts,cex.lab=(2),cex.axis=(1.5))
    usr <- par("usr")
    text(usr[2], usr[4], grain, adj = c(1.2, 1.2), cex = 1.2)
    if (save==TRUE){
      filename=file.path(out_dir, paste0( grain,"_zoom", ".pdf"))
      dev.copy2pdf(file = filename, width = 8, height = 6)
      }

}

make_grain_wetherill <- function(df, grain, prop, pts, save) {
  
  grain_data <- splt_data(df, grain, prop)
  
  if (save) {
    filename <- file.path(out_dir, paste0(grain, "full.pdf"))
    pdf(filename, width = 8, height = 6)
  }
  
  concordia(
    grain_data, type = 1, anchor = c(2,0),
    ellipse.fill = alpha('#de2d26', 0.6),
    ticks = 10, concordia.col = 'black',
    show.numbers = FALSE, show.age = 2,
    oerr = 2, common.Pb = 0, omit = pts,
    cex.lab = 2, cex.axis = 1.5
  )
  
  usr <- par("usr")
  text(usr[2], usr[4], grain, adj = c(1.2, 1.2), cex = 1.2)
  
  if (save) dev.off()
}


# calculate the 7/6 age for acropolis

age_76<- function(df,sample,prop){
  grain_data <- splt_data(df, sample, prop)
  indiv_age <- age(grain_data, method = "Pb207-Pb206", oerr = 2, type = 3, common.Pb = 0, isochron = TRUE)

  age_val <- indiv_age$par[1]
  pb_int<-indiv_age$par[2]
  mswd <- indiv_age$mswd
  sigma <- indiv_age$err[1]
  ci <- indiv_age$err[2]
  
  ages <- data.frame(
    sample = sample,
    age = age_val,
    mswd = mswd,
    sigma = sigma,
    ci = ci)
  
  return(ages)
}

# This is to plot limited regions of the TW concordias for the paper 

plot_inset<- function(df, grain,prop,xlim,ylim, hide,pts,save){
  
  
  if(hide==FALSE){
    grain_data<-splt_data(df,grain,prop)
    
    concordia(grain_data,type=2,ellipse.fill = alpha('#de2d26', 0.6), ticks = 10,
              show.numbers = FALSE,show.age = 2, oerr = 2, common.Pb = 0,cex.lab=(2),cex.axis=(1.5),xlim=xlim,
              ylim=ylim)
    usr <- par("usr")
    text(usr[2], usr[4], grain, adj = c(1, 1), cex = 1.2)
    
  }
  
  if(hide==TRUE){
    grain_data<-splt_data(df,grain,prop)
    print(rownames(grain_data))
    concordia(grain_data,type=2,ellipse.fill = alpha('#de2d26', 0.6), ticks = 10,
              show.numbers = FALSE,show.age = 2, oerr = 2, common.Pb = 0,omit=pts,cex.lab=(2),cex.axis=(1.5),xlim=xlim,
              ylim=ylim)
    
    usr <- par("usr")
    text(usr[2], usr[4], grain, adj = c(1, 1), cex = 1.2)
    if (save==TRUE){
      filename=file.path(out_dir, paste0( grain,"inset", ".pdf"))
      dev.copy2pdf(file = filename, width = 8, height = 6)
    }
    
    
  }
  
}

## This is for filtering based upon the assigned mineral name (titano vs hematite )
filter_composition<- function(upb,comp){

  filt_df<-upb%>% filter(sample_grp %in% comp$Grain)
  print(unique(filt_df$sample_grp))
  return(filt_df)
}
  







####### Main Body of code #########



out_dir <- "~/GitHub/Sed_Hematite_UPb/figures/geochron/single_grain"


#### Read in the data ##### 
upb_data <- read.csv("data/geochron/plotting/upb_data/2025_day1_fin.csv")
upb_data2 <- read.csv("data/geochron/plotting/upb_data/2025_day2_fin.csv")
upb_data3<-read.csv("data/geochron/plotting/upb_data/2025_day3_fin.csv")
upb_data4<-read.csv("data/geochron/plotting/upb_data/2025_day4_fin.csv")
upb_data5<-read.csv("data/geochron/plotting/upb_data/2025_day5.csv")
upb_data_old <- read.csv("data/geochron/plotting/upb_data/2024_dh.csv")
upb_data_old_02 <- read.csv("data/geochron/plotting/upb_data/2024_dh_02.csv")


# Compositions 

hem_grains<-read.csv("data/eds/class/total_hematite_v2.csv")
titano_grains<-read.csv("data/eds/class/total_titano_hematite_v2.csv")
hematite_grains<-unique(rbind(hem_grains,titano_grains))



##### Filter by Pb spread #######

# Set the Pb ratio filter here - 

sample_spread=0.04

### Filter by user specified filter 
day_01<-filter_by_spread(upb_data,name_filter='LRa',thrs=sample_spread,
                         type='Pb',scale='absolute')
day_02<-filter_by_spread(upb_data2,name_filter='LRa',thrs=sample_spread,
                         type='Pb',scale='absolute')

day_03<-filter_by_spread(upb_data3,name_filter='BRe',thrs=sample_spread,
                         type='Pb',scale='absolute')

day_04<-filter_by_spread(upb_data4,name_filter='BRc',thrs=sample_spread,
                         type='Pb',scale='absolute')

day_05<-filter_by_spread(upb_data5,name_filter='BRc',thrs=sample_spread,
                         type='Pb',scale='absolute')

apr_24<-filter_by_spread(upb_data_old,name_filter='BRc',thrs=sample_spread,
                         type='Pb',scale='absolute')

apr_24_02<-filter_by_spread(upb_data_old_02,name_filter='BRc',thrs=sample_spread,
                            type='Pb',scale='absolute')



##### Acropolis #####
acropolis_01<-sample_spread_filter(upb_data,0.0,'Acropolis')
acropolis_02<-sample_spread_filter(upb_data2,0.0,'Acropolis')
acropolis_03<-sample_spread_filter(upb_data3,0.0,'Acropolis')
acropolis_04<-sample_spread_filter(upb_data4,0.0,'Acropolis')
acropolis_05<-sample_spread_filter(upb_data5,0.0,'Acropolis')



#### Filter by mineral####
day_01_hem<-filter_composition(day_01,hem_grains)
day_01_titano<-filter_composition(day_01,titano_grains)

day_02_hem<-filter_composition(day_02,hem_grains)
day_02_titano<-filter_composition(day_02,titano_grains)

day_03_hem<-filter_composition(day_03,hem_grains)
day_03_titano<-filter_composition(day_03,titano_grains)

day_04_hem<-filter_composition(day_04,hem_grains)
day_04_titano<-filter_composition(day_04,titano_grains)


day_05_hem<-filter_composition(day_05,hem_grains)
day_05_titano<-filter_composition(day_05,titano_grains)

apr_24_hem<-filter_composition(apr_24,hem_grains)
apr_24_titano<-filter_composition(apr_24, titano_grains)

apr_24_02_hem<-filter_composition(apr_24_02,hem_grains)
apr_24_02_titano<-filter_composition(apr_24_02, titano_grains)

##### Plot sample groups from spread filter to triage for additional interpretation ####
#### This region of the code is just to get an overview by day- all the plots of grains passing the 
### filtering threshold will plot all at once so if you don't want that skip

prop_val <- FALSE

# First hematite 
plot_by_sample(apr_24_hem,prop_val,save=FALSE)
plot_by_sample(apr_24_02_hem,prop_val,save=FALSE)
plot_by_sample(day_01_hem,prop_val,save=FALSE)
plot_by_sample(day_02_hem,prop_val,save=FALSE)
plot_by_sample(day_03_hem,prop_val,save=FALSE)
plot_by_sample(day_04_hem,prop_val,save=FALSE)
plot_by_sample(day_05_hem,prop_val,save=FALSE)

plot_by_sample(day_01_test,prop_val,save=FALSE)

# Second the titanohematite
plot_by_sample(apr_24_titano,prop_val,save=FALSE)
plot_by_sample(apr_24_02_titano,prop_val,save=FALSE)
plot_by_sample(day_01_titano,prop_val,save=FALSE)

plot_by_sample(day_02_titano,prop_val,save=FALSE)
plot_by_sample(day_03_titano,prop_val,save=FALSE)
plot_by_sample(day_04_titano,prop_val,save=FALSE)
plot_by_sample(day_05_titano,prop_val,save=FALSE)

plot_wetherill_sample(acropolis_01,prop_val,save=FALSE)
plot_wetherill_sample(acropolis_02,prop_val,save=FALSE)
plot_wetherill_sample(acropolis_03,prop_val,save=FALSE)
plot_wetherill_sample(acropolis_04,prop_val,save=FALSE)

plot_wetherill_sample(day_02,prop_val,save=FALSE)








##### Inspect the samples needing additional pruning ####

plot_grains(apr_24,'BRc_grain_33',prop_val,hide=TRUE,pts=c())
plot_grains(apr_24_02,'BRc_grain_93',prop_val,hide=TRUE,pts=c(3))
plot_grains(day_02,'LRa02_grain_56',prop_val,hide=TRUE,pts=c())
plot_grains(day_03,'BRe_grain_05',prop_val,hide=TRUE,pts=c(3,4,5))
plot_grains(day_03,'BRe_grain_101',prop_val,hide=TRUE,pts=c(4))
plot_grains(day_03,'BRe_grain_105',prop_val,hide=TRUE,pts=c(),zoom=TRUE)
plot_grains(day_03,'BRe_grain_100',prop_val,hide=TRUE,pts=c(),zoom=TRUE)
plot_grains(day_04,'BRc02_grain_06',prop_val,hide=TRUE,pts=c(2),zoom=TRUE)

plot_grains(day_01,'LRa01_grain_48',prop_val,hide=TRUE,pts=c())

plot_grains(day_02,'LRa02_grain_50',prop_val,hide=TRUE,pts=c())

plot_grains(day_02,'LRa02_grain_38',prop_val,hide=TRUE,pts=c())



plot_grains(day_02,'LRa02_grain_94',prop_val,hide=TRUE,pts=c())




plot_grains(apr_24_titano,'BRc_grain_22',prop_val,hide=TRUE,pts=c(1))
plot_grains(apr_24_02_titano,'BRc_grain_93',prop_val,hide=TRUE,pts=c(1,3,7))            
plot_grains(apr_24_02_titano,'BRc_grain_51',prop_val,hide=TRUE,pts=c(6,7,10))
plot_grains(day_01_titano,'LRa01_grain_21',prop_val,hide=TRUE,pts=c(1))
plot_grains(day_01_titano,'LRa01_grain_12',prop_val,hide=TRUE,pts=c(1,2,3,4))
plot_grains(day_02_titano,'LRa02_grain_90',prop_val,hide=TRUE,pts=c(5))
plot_grains(day_02_titano,'LRa02_grain_70',prop_val,hide=TRUE,pts=c(6))
plot_grains(day_04,'BRc02_grain_97',prop_val,hide=TRUE,pts=c(1,4,5,14,15),zoom=TRUE) 
plot_grains(day_04,'BRc02_grain_06',prop_val,hide=TRUE,pts=c(2),zoom=TRUE)

 

plot_grains_wetherill(acropolis_01,'Acropolis',prop = prop_val,hide = TRUE,pts=c(12,14,15,16,19),zoom=TRUE)
plot_grains_wetherill(acropolis_03,'Acropolis',prop = prop_val,hide = TRUE,pts=c(19,2,22,15,9),zoom=TRUE)
plot_grains_wetherill(acropolis_04,'Acropolis',prop = prop_val,hide = TRUE,pts=c(),zoom=TRUE)


lra_grain<-splt_data(day_02_titano,'LRa02_grain_90','FALSE')
lra_age<-age(lra_grain,method = "U238-Pb206", oerr = 3, 
             type = 3, common.Pb = 2, isochron = FALSE)
lra_age$err
# Filter out bad spots and save
make_grain_plots(day_04,'BRc02_grain_97',prop_val, pts=c(1,4,5),save=TRUE)
make_grain_plots(day_01,'LRa01_grain_12',prop_val, pts=c(1,2,3,4),save=TRUE)
make_grain_plots(day_02,'LRa02_grain_56',prop_val, pts=c(),save=TRUE)
make_grain_plots(day_01,'LRa01_grain_48',prop_val, pts=c(),save=TRUE)

make_grain_plots(day_03,'BRe_grain_05',prop_val, pts=c(),save=TRUE)
make_grain_plots(day_03,'BRe_grain_105',prop_val, pts=c(),save=TRUE)
make_grain_plots(day_04,'BRc02_grain_35',prop_val, pts=c(),save=TRUE)
make_grain_plots(day_04,'BRc02_grain_06',prop_val, pts=c(2),save=TRUE)
make_grain_plots(apr_24_02_titano,'BRc_grain_51',prop_val,hide=TRUE,pts=c(6,7,10))
make_grain_plots(day_05,'BRc_grain_109',prop_val, pts=c(8,10),save=TRUE)



make_grain_plots(day_04,'BRc02_grain_61',prop_val, pts=c(),save=TRUE)


make_grain_wetherill(acropolis_01,'Acropolis',prop_val,pts=c(),save=TRUE)
make_grain_wetherill(acropolis_03,'Acropolis',prop_val,pts=c(,save=TRUE)
make_grain_wetherill(acropolis_04,'Acropolis',prop_val,pts=c(),save=TRUE)

#### Plot the insets if need to change aspect ratio ####
plot_inset(day_01,'LRa01_grain_48',xlim=c(0,0.7),ylim=c(0.79,0.86),prop_val,hide=TRUE,pts=c(),save=TRUE)
plot_inset(day_04,'BRc02_grain_06',xlim=c(0.2,0.8),ylim=c(0.76,0.82),prop_val,hide=TRUE,pts=c(2),save=TRUE)
plot_inset(day_01,'LRa01_grain_12',xlim=c(0.9,1.9),ylim=c(0.62,0.72),prop_val,hide=TRUE,pts=c(1,2,3,4),save=TRUE)
plot_inset(day_04,'BRc02_grain_97',xlim=c(1,3),ylim=c(0.5,0.7),prop_val,hide=TRUE,pts=c(1,4,5),save=TRUE)
plot_inset(day_05,'BRc_grain_137',prop_val,xlim=c(0.7,2.3),ylim=c(0.6,0.76),hide=TRUE,pts=c(3,6,10,12,14))


# Inset for bad grains 
plot_inset(day_05,'BRc_grain_111',xlim=c(0.2,0.8),ylim=c(0.76,0.82),prop_val,hide=TRUE,pts=c(),save=TRUE)
# Inset for bad grains 
plot_inset(day_04,'BRc02_grain_61',xlim=c(0.2,0.8),ylim=c(0.79,0.85),prop_val,hide=TRUE,pts=c(),save=TRUE)


#### recombine for export####
apr_24<-unique(rbind(apr_24_hem,apr_24_titano))
apr_24_02<-unique(rbind(apr_24_02_hem,apr_24_02_titano))
day_01<-unique(rbind(day_01_hem,day_01_titano))
day_02<-unique(rbind(day_02_hem,day_02_titano))
day_03<-unique(rbind(day_03_hem,day_03_titano))
day_04<-unique(rbind(day_04_hem,day_04_titano))
day_05<-unique(rbind(day_05_hem,day_05_titano))



##### Filter out data by importing csvs for each day ####
apr_24_pts<-read.csv('data/geochron/plotting/upb_data/omit/0424_01_omit.csv',stringsAsFactors = FALSE)
apr_24<-remove_pts(apr_24,apr_24_pts)

apr_24_02_pts<-read.csv('data/plotting/geochron/plotting/upb_data/omit/0424_02_omit.csv',stringsAsFactors = FALSE)
apr_24_02<-remove_pts(apr_24_02,apr_24_02_pts)

day_01_pts<-read.csv('data/plotting/geochron/plotting/upb_data/omit/day01_omit.csv',stringsAsFactors = FALSE)
day_01<-remove_pts(day_01,day_01_pts)

day_02_pts<-read.csv('data/plotting/geochron/plotting/upb_data/omit/day02_omit.csv',stringsAsFactors = FALSE)
day_02<-remove_pts(day_02,day_02_pts)


day_03_pts<-read.csv('data/plotting/geochron/plotting/upb_data/omit/day03_omit.csv',stringsAsFactors = FALSE)
day_03<-remove_pts(day_03,day_03_pts)

day_04_pts<-read.csv('data/plotting/geochron/plotting/upb_data/omit/day04_omit.csv')
day_04<-remove_pts(day_04,day_04_pts)

day_05_pts<-read.csv('data/geochron/plotting/upb_data/omit/day05_omit.csv')
day_05<-remove_pts(day_05,day_05_pts)






