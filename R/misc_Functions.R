# Generate survival plots with p-values and median OS.
plotSurvival <- function(fit, ylim, data, palette = 'jco', hr = NULL){
    
    # Generate survival plot.
    x <- survminer::ggsurvplot(
        fit = fit,
        pval = F,
        size = .825,
        break.time.by = 10,
        break.y.by = .2,
        palette = palette,
        risk.table = T,
        tables.height = .3,
        xlab = 'Time (in months)',
        axes.offset = F,
        ylim = c(0, 1.05),
        xlim = c(0, ylim), strata.labels = 'c',
        risk.table.col = 'strata', censor.shape = '+',
        fontsize = 3,
        risk.table.title = 'No. at risk',
        ggtheme = ggplot2::theme(
            legend.position = 'bottom',
            legend.direction = 'horizontal',
            text = ggplot2::element_text(size=9, family='Arial', face = 'bold'),
            axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
            axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
            panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
            panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
            legend.text = ggtext::element_markdown()
        )
    )
    
    # Add the log-rank p-value.
    p.logrank <- survminer::surv_pvalue(fit = fit, method = 'log-rank', data = data, test.for.trend = F)
    x$plot <- x$plot + ggplot2::annotate('text', x = max(x$data.survplot$time) * .75, y = 1, label = paste0('log-rank: ', p.logrank$pval.txt), size = 2.5)
    
    # Add HR (if two groups)
    if(!is.null(hr)){
        
        HR.CI <- round(summary(hr)$conf.int, 2)
        HR.p <- round(summary(hr)$waldtest[[3]], 2)
        HR.CI <- sprintf('HR (95%% CI): %s (%s - %s)', HR.CI[[1]], HR.CI[[3]], HR.CI[[4]])
        x$plot <- x$plot + ggplot2::annotate('text', x = max(x$data.survplot$time) * .75, y = .9, label = HR.CI, size = 2.5)
    }
    
    # Add the median OS.
    medianOS <- x$data.survplot %>%
        dplyr::group_by(strata) %>%
        dplyr::summarise(
            medianOS = round(median(time, na.rm = T), 2),
            label = sprintf('%s - %s mo.', unique(strata), medianOS)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(-medianOS)
    
    x$plot <- x$plot + ggplot2::annotate('text', x = max(x$data.survplot$time) * .75, y = .75, label = paste0('Median OS (Desc.):\n', paste(medianOS$label, collapse = '\n')), size = 2.5)
    
    # Remove legends.
    x$plot <- x$plot + theme_Job + theme(legend.position = 'none')
    x$table <- x$table + theme_Job + theme(legend.position = 'none')
    
    # Add perc. scaling.
    x$plot <- x$plot + ggplot2::scale_y_continuous(labels = scales::percent_format())
    
    return(x)
}

plotHR <- function(data, withQ = F){
    x <- data %>% 
        gtsummary::tbl_regression(
            exponentiate = T, 
            add_estimate_to_reference_rows = T,
        ) %>%
        gtsummary::add_n() %>% 
        gtsummary::add_global_p() %>% 
        gtsummary::add_nevent() %>% 
        gtsummary::bold_p() %>% 
        gtsummary::bold_labels() %>% 
        gtsummary::italicize_levels() %>% 
        gtsummary::sort_p() %>% 
        bstfun::add_inline_forest_plot(header = '', spec_pointrange.args = list(lim = c(-3, 3), width = 550, cex = 1, col = 'black', pch = 1))
    
    if(withQ){
        x <- x %>% gtsummary::add_q()
    }
    
    x %>% bstfun::as_ggplot()
    
}
