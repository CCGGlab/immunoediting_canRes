table_apply_fun <- function(tab, row=NULL, column=NULL, fun)
{
  tabGrob <- ggpubr:::get_tablegrob(tab)
  
  # Code added to support easily applying the function on all rows and / columns
  allcols <- tabGrob$layout %>% filter(name == "core-fg") %>% pull(l) %>% unique
  allrows <- tabGrob$layout %>% filter(name == "core-fg") %>% pull(t) %>% unique
  
  if (is.null(column))
    column <- allcols
  if (is.null(row))
    row <- allrows
  
  cells <- expand.grid(row = row, column = column)
  for(i in 1:nrow(cells)){
    tc <- ggpubr:::.find_cell(tabGrob, cells$row[i], cells$column[i], "core-fg")
    tabGrob$grobs[tc][[1]] <- fun(tabGrob$grobs[tc][[1]])
  }
  
  ggpubr:::tab_return_same_class_as_input(tabGrob, input = tab)
}
