#Accuracy
top_genes <- importance_df_accuracy |> dplyr::arrange(desc(MeanDecreaseAccuracy)) |> dplyr::pull(Gene) |> head(100)
cor_mat <- cor(train_data_clean[, top_genes[top_genes%in% colnames(train_data_clean)]])
heatmap(cor_mat)

#Gini
top_genes <- importance_df_gini |> dplyr::arrange(desc(MeanDecreaseGini)) |> dplyr::pull(Gene) |> head(100)
cor_mat <- cor(train_data_clean[, top_genes[top_genes%in% colnames(train_data_clean)]])
heatmap(cor_mat)