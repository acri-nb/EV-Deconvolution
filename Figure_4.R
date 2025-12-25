library(fmsb)
library(dplyr)

# Ground truth data
ground_truth_data <- data.frame(
  Mix = c("Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6"),
  Breast = c(0.15, 0, 0.15, 0, 0.15, 0),
  Blymphocyte = c(0.05, 0, 0.05, 0, 0.05, 0.1),
  `Blymphocyte.stim` = c(0.1, 0, 0.1, 0, 0.1, 0.2),
  Kidney = c(0.5, 0, 0.5, 0.1, 0.5, 0),
  Colon = c(0.2, 0, 0.2, 0.05, 0.2, 0),
  Myeloid = c(0, 0.1, 0, 0, 0, 0.2),
  Lung = c(0, 0.05, 0, 0.15, 0, 0),
  Pancreas = c(0, 0.15, 0, 0.5, 0, 0),
  Platelets = c(0, 0.2, 0, 0, 0, 0.2),
  Prostate = c(0, 0.5, 0, 0.2, 0, 0),
  TLymphocyte = c(0, 0, 0, 0, 0, 0.1),
  `Tlymphocyte.Stim` = c(0.1, 0, 0.1, 0, 0.1, 0.2)
)

# Convert to proportions
row_sums <- rowSums(ground_truth_data[, -1])
ground_truth_proportions <- ground_truth_data
for (i in 2:ncol(ground_truth_proportions)) {
  ground_truth_proportions[, i] <- ground_truth_proportions[, i] / row_sums
}

# Method data (CIBERSORT, MuSiC, etc.)
cibersort_data <- data.frame(
  Mix = c("Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6"),
  Blymphocyte = c(0.03544808, 0, 0.03252584, 0, 0.03250260, 0.08171565),
  Blymphocyte.stim = c(0.08664369, 0.00126137, 0.08088334, 0, 0.08271591, 0.09938276),
  Breast = c(0.11235536, 0.01629801, 0.09664910, 0, 0.09964531, 0),
  Colon = c(0.25059464, 0, 0.23872318, 0.02845874, 0.23804884, 0.06797212),
  Kidney = c(0.43977047, 0.00179195, 0.45368759, 0.09968292, 0.45652048, 0.00471870),
  Lung = c(0, 0, 0, 0.09055198, 0, 0),
  Myeloid = c(0.00849479, 0.11615454, 0.01207317, 0.03619246, 0.01252438, 0.18837339),
  Pancreas = c(0, 0.09053287, 0, 0.35302824, 0, 0.07521662),
  Platelets = c(0.00341864, 0.26163745, 0.00324551, 0.00125320, 0.00341519, 0.22299663),
  Prostate = c(0.05182007, 0.44337890, 0.05666456, 0.33556718, 0.05160016, 0),
  TLymphocyte = c(0.01145425, 0.06894490, 0.02554770, 0.05526527, 0.02302712, 0.16715631),
  Tlymphocyte.Stim = c(0, 0, 0, 0, 0, 0.09246781)
)

music_data <- data.frame(
  Mix = c("Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6"),
  Blymphocyte = c(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000),
  Blymphocyte.stim = c(0.1521, 0.0292, 0.1518, 0.0000, 0.1518, 0.3735),
  Breast = c(0.1379, 0.0895, 0.1409, 0.1444, 0.1455, 0.0000),
  Colon = c(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000),
  Kidney = c(0.4850, 0.0000, 0.4841, 0.1926, 0.4845, 0.0118),
  Lung = c(0.0386, 0.0000, 0.0405, 0.0272, 0.0388, 0.0000),
  Myeloid = c(0.0000, 0.0334, 0.0000, 0.0000, 0.0000, 0.0840),
  Pancreas = c(0.1832, 0.6284, 0.1785, 0.4941, 0.1753, 0.0000),
  Platelets = c(0.0031, 0.2194, 0.0042, 0.0000, 0.0041, 0.3105),
  Prostate = c(0.0000, 0.0000, 0.0000, 0.1417, 0.0000, 0.0000),
  TLymphocyte = c(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000),
  Tlymphocyte.Stim = c(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.2202)
)

dual_simplex_data <- data.frame(
  Mix = c("Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6"),
  Blymphocyte = c(2.976693e-18, 7.238377e-18, 3.810873e-18, -7.684857e-17, 1.201748e-18, 3.167624e-01),
  Blymphocyte.stim = c(7.574810e-17, -1.390507e-17, 1.199754e-16, -6.873815e-17, 7.879790e-17, -4.377621e-17),
  Breast = c(2.211915e-01, -1.902651e-18, 2.015607e-01, 3.899711e-02, 2.230418e-01, 8.847443e-18),
  Colon = c(1.738153e-17, 5.949081e-17, -3.138531e-17, 2.658724e-17, -2.370326e-17, 3.372543e-17),
  Kidney = c(7.788085e-01, -9.151554e-20, 7.789732e-01, 5.110100e-01, 7.769582e-01, 2.137894e-02),
  Lung = c(2.292003e-16, -2.691641e-18, 1.467332e-16, -1.953111e-16, 1.623951e-16, -3.679393e-18),
  Myeloid = c(-3.110588e-17, -1.203341e-17, 3.859512e-18, -9.090631e-18, 3.661281e-18, 1.622476e-01),
  Pancreas = c(9.840948e-20, 2.812259e-01, 1.946611e-02, 4.457316e-01, -8.456681e-20, 9.897054e-18),
  Platelets = c(1.965248e-19, 4.597098e-01, -2.979761e-18, -2.174889e-18, -7.746821e-18, 4.996111e-01),
  Prostate = c(-4.045234e-16, 2.590643e-01, -3.853790e-16, 4.261314e-03, -3.972425e-16, -4.676349e-17),
  TLymphocyte = c(1.846640e-17, 8.554129e-17, -2.019911e-17, -6.120612e-17, -9.974101e-18, 3.153791e-18),
  Tlymphocyte.Stim = c(4.555504e-18, -1.226786e-16, -1.787359e-16, -7.578306e-18, -1.620348e-16, -3.279305e-17)
)

nmf_data <- data.frame(
  Mix = c("Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6"),
  Kidney = c(0.457979119488629, 0.185977102488706, 0.531021251206337, 0.502372808205509, 0.533846250843462, 0),
  Prostate = c(0.148089221818423, 0.143953065170294, 0.128478892224605, 0.369368737367811, 0.128807232839981, 0),
  Blymphocyte = c(0.291652417773992, 0, 0.279720122165664, 0.0901060282196631, 0.281617015935284, 0.642201452023517),
  Lung = c(0.0850678482028392, 0, 0.0400326835047353, 0, 0.033927058695742, 0),
  Platelets = c(0.0172113927161165, 0.670069832341001, 0.0207470508986588, 0.0381524262070174, 0.0218024416855307, 0.357798547976483),
  Blymphocyte.stim = c(0, 0, 0, 0, 0, 0),
  Breast = c(0, 0, 0, 0, 0, 0),
  Colon = c(0, 0, 0, 0, 0, 0),
  Myeloid = c(0, 0, 0, 0, 0, 0),
  Pancreas = c(0, 0, 0, 0, 0, 0),
  TLymphocyte = c(0, 0, 0, 0, 0, 0),
  Tlymphocyte.Stim = c(0, 0, 0, 0, 0, 0)
)

nlls_data <- data.frame(
  Mix = c("Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6"),
  Blymphocyte = c(0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.01346631),
  Blymphocyte.stim = c(0.14497517, 0.06518118, 0.14463538, 0.07297522, 0.14496675, 0.17008175),
  Breast = c(0.20239329, 0.07195102, 0.20496013, 0.06024032, 0.20700451, 0.04195438),
  Colon = c(0, 0, 0, 0, 0, 0),
  Kidney = c(0.20581673, 0.00000000, 0.20448794, 0.06222798, 0.20492304, 0.00000000),
  Lung = c(0.114395529, 0.075318773, 0.115758372, 0.173507490, 0.114689313, 0.002686254),
  Myeloid = c(0.00000000, 0.07400891, 0.00000000, 0.00000000, 0.00000000, 0.14655733),
  Pancreas = c(0.2176951, 0.2787699, 0.2137560, 0.4543043, 0.2124891, 0.2382445),
  Platelets = c(0.02487181, 0.21880772, 0.02554621, 0.01358061, 0.02555993, 0.21202936),
  Prostate = c(0.0000000, 0.1941958, 0.0000000, 0.1175246, 0.0000000, 0.0000000),
  TLymphocyte = c(0, 0, 0, 0, 0, 0),
  Tlymphocyte.Stim = c(0.08985239, 0.02176675, 0.09085595, 0.04563950, 0.09036731, 0.17498015)
)

nlls_quadprog_data <- data.frame(
  Mix = c("Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6"),
  Blymphocyte = c(0.13890662, 0.01790041, 0.13438085, 0.00000000, 0.13809594, 0.15230924),
  Blymphocyte.stim = c(0.01885689, 0.11195282, 0.02396899, 0.09890843, 0.01059959, 0.11621152),
  Breast = c(0.15009062, 0.05285325, 0.12813621, 0.06800441, 0.13039702, 0.06593136),
  Colon = c(0.18793321, 0.03400344, 0.18371456, 0.06060474, 0.18506425, 0.03504738),
  Kidney = c(0.23002136, 0.00000000, 0.22326926, 0.08466211, 0.23301971, 0.00000000),
  Lung = c(0.06017492, 0.08174585, 0.07365307, 0.15710046, 0.06941715, 0.04274868),
  Myeloid = c(0.05374678, 0.10938500, 0.05719799, 0.02843053, 0.05552310, 0.13377374),
  Pancreas = c(0.00000000, 0.09745456, 0.02085993, 0.18935251, 0.01879435, 0.06185263),
  Platelets = c(0.03130424, 0.13532285, 0.02645816, 0.04282437, 0.02304654, 0.12940003),
  Prostate = c(0.000000000, 0.231680540, 0.000000000, 0.120704113, 0.000000000, 0.006297306),
  TLymphocyte = c(0.11643794, 0.03882258, 0.11013243, 0.03525185, 0.10816389, 0.13252163),
  Tlymphocyte.Stim = c(0.01252742, 0.08887870, 0.01822855, 0.11415648, 0.02787847, 0.12390649)
)

hspe_data <- data.frame(
  Mix = c("Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6"),
  Blymphocyte = c(0.0615462915975005, 2.10986670974832e-14, 0.056581668354361, 1.63333935287624e-14, 0.0538422302761945, 0.109253055838156),
  Blymphocyte.stim = c(0.17205189564674, 8.9529327722133e-15, 0.149085653208496, 2.12357239506118e-14, 0.158076349966032, 0.565480232581023),
  Breast = c(0.159174908896492, 1.12054921216177e-12, 0.157809904032517, 2.37654512353938e-12, 0.168462474830569, 8.15028066494612e-14),
  Colon = c(0.184944097814528, 0.000860432365042374, 0.215606708016647, 0.0440677338966095, 0.216085511553033, 0.0017857447927055),
  Kidney = c(0.374188781259241, 0.000276500335019246, 0.369621003708629, 0.0638667986895588, 0.359589944467675, 0.000600596018706707),
  Lung = c(0.00113345682398526, 0.0555962521047241, 0.00253299657465361, 0.0691174535927729, 0.00392343153901516, 0.0318630911936676),
  Myeloid = c(0.000410236518530606, 0.0405496285468119, 0.000551736347792372, 0.000350613878038412, 0.000542608602068294, 0.0826691200391667),
  Pancreas = c(0.0419007094640993, 0.125321093618446, 0.0445982183455745, 0.480411225036451, 0.0280126639522533, 0.000201941790806037),
  Platelets = c(8.29940111804725e-06, 0.0500390217826638, 1.40926938417996e-05, 2.48108732758093e-06, 1.23863298025725e-05, 0.0534074739655798),
  Prostate = c(2.25871144600948e-11, 0.726862588963815, 0.0016033206481104, 0.341233007697635, 0.00695016108583158, 1.19323243989806e-12),
  TLymphocyte = c(0.0046413225549792, 2.01007047484897e-05, 0.00199469806852641, 0.0003042093938911, 0.00210970378859576, 0.0726900406644244),
  Tlymphocyte.Stim = c(1.98289893473012e-13, 0.000474381577578971, 8.49851697974119e-13, 0.000646476725300891, 0.00239253360892957, 0.08204870311449)
)

# Calculate metrics function
calculate_metrics <- function(predicted_data, ground_truth_data, method_name) {
  pred_matrix <- as.matrix(predicted_data[, -1])
  truth_matrix <- as.matrix(ground_truth_data[, -1])
  
  pred_matrix <- apply(pred_matrix, 2, as.numeric)
  truth_matrix <- apply(truth_matrix, 2, as.numeric)
  
  rmse <- sqrt(mean((pred_matrix - truth_matrix)^2))
  mae <- mean(abs(pred_matrix - truth_matrix))
  correlation <- cor(as.vector(pred_matrix), as.vector(truth_matrix))
  
  ss_res <- sum((as.vector(pred_matrix) - as.vector(truth_matrix))^2)
  ss_tot <- sum((as.vector(truth_matrix) - mean(as.vector(truth_matrix)))^2)
  r_squared <- 1 - (ss_res / ss_tot)
  
  threshold <- 0.01
  pred_binary <- pred_matrix > threshold
  truth_binary <- truth_matrix > threshold
  
  tp <- sum(pred_binary & truth_binary)
  fp <- sum(pred_binary & !truth_binary)
  tn <- sum(!pred_binary & !truth_binary)
  fn <- sum(!pred_binary & truth_binary)
  
  sensitivity <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  specificity <- ifelse(tn + fp > 0, tn / (tn + fp), 0)
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  f1_score <- ifelse(precision + sensitivity > 0, 2 * (precision * sensitivity) / (precision + sensitivity), 0)
  
  cosine_sim <- mean(apply(cbind(pred_matrix, truth_matrix), 1, function(row_data) {
    pred_row <- row_data[1:ncol(pred_matrix)]
    truth_row <- row_data[(ncol(pred_matrix)+1):(2*ncol(pred_matrix))]
    sum(pred_row * truth_row) / (sqrt(sum(pred_row^2)) * sqrt(sum(truth_row^2)))
  }))
  
  return(data.frame(
    Method = method_name,
    RMSE = rmse,
    MAE = mae,
    Correlation = correlation,
    R_squared = r_squared,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1_Score = f1_score,
    Cosine_Similarity = cosine_sim
  ))
}

# Calculate metrics for all methods
all_metrics <- bind_rows(
  calculate_metrics(cibersort_data, ground_truth_proportions, "CIBERSORT"),
  calculate_metrics(music_data, ground_truth_proportions, "MuSiC"),
  calculate_metrics(dual_simplex_data, ground_truth_proportions, "Dual Simplex"),
  calculate_metrics(nmf_data, ground_truth_proportions, "NMF"),
  calculate_metrics(nlls_data, ground_truth_proportions, "NLLS"),
  calculate_metrics(nlls_quadprog_data, ground_truth_proportions, "NLLS + Quadprog"),
  calculate_metrics(hspe_data, ground_truth_proportions, "HSPE")
)

# Prepare radar data
radar_data <- all_metrics
radar_data$RMSE_Accuracy <- 1 - (radar_data$RMSE / max(radar_data$RMSE))
radar_data$MAE_Accuracy <- 1 - (radar_data$MAE / max(radar_data$MAE))

radar_data <- radar_data[, c("Method", "RMSE_Accuracy", "MAE_Accuracy", "Correlation", 
                             "R_squared", "Sensitivity", "Specificity", "Precision", 
                             "F1_Score", "Cosine_Similarity")]

# Prepare for fmsb
radar_matrix <- as.matrix(radar_data[, -1])
rownames(radar_matrix) <- radar_data$Method

radar_for_plot <- rbind(
  rep(1, ncol(radar_matrix)),
  rep(0, ncol(radar_matrix)),
  radar_matrix
)
radar_for_plot <- as.data.frame(radar_for_plot)

# Create radar plot
colors <- rainbow(nrow(radar_data))

pdf("deconvolution_radar_simple.pdf", width = 12, height = 8)
par(mfrow = c(3, 3), mar = c(1, 1, 2, 1))

# Individual method plots
for (i in 1:nrow(radar_data)) {
  single_method_data <- as.data.frame(radar_for_plot[c(1, 2, i + 2), ])
  radarchart(single_method_data,
             axistype = 1,
             pcol = colors[i],
             pfcol = scales::alpha(colors[i], 0.3),
             plwd = 3,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "black",
             caxislabels = seq(0, 1, 0.2),
             cglwd = 1.2,
             vlcex = 0.9,
             title = radar_data$Method[i])
}

# All methods comparison overlay
radarchart(as.data.frame(radar_for_plot),
           axistype = 1,
           pcol = colors,
           pfcol = scales::alpha(colors, 0.1),
           plwd = 3,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "black",
           caxislabels = seq(0, 1, 0.2),
           cglwd = 1.2,
           vlcex = 0.9,
           title = "All Methods Comparison")

# Legend panel
plot.new()
legend("center",
       legend = radar_data$Method,
       col = colors,
       lty = 1,
       lwd = 3,
       cex = 1.0,
       bty = "n")

dev.off()

# Also save as PNG at 600 DPI
png("deconvolution_radar_simple.png", width = 3600, height = 2700, res = 600)
par(mfrow = c(3, 3), mar = c(1, 1, 2, 1), family = "sans", font = 2, font.lab = 2, font.axis = 2)

for (i in 1:nrow(radar_data)) {
  single_method_data <- as.data.frame(radar_for_plot[c(1, 2, i + 2), ])
  radarchart(single_method_data,
             axistype = 1,
             pcol = colors[i],
             pfcol = scales::alpha(colors[i], 0.3),
             plwd = 3,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "black",
             caxislabels = seq(0, 1, 0.2),
             cglwd = 1.2,
             vlcex = 0.9,
             title = radar_data$Method[i])
}

# All methods comparison overlay
radarchart(as.data.frame(radar_for_plot),
           axistype = 1,
           pcol = colors,
           pfcol = scales::alpha(colors, 0.1),
           plwd = 3,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "black",
           caxislabels = seq(0, 1, 0.2),
           cglwd = 1.2,
           vlcex = 0.9,
           title = "All Methods Comparison")

plot.new()
legend("center",
       legend = radar_data$Method,
       col = colors,
       lty = 1,
       lwd = 3,
       cex = 1.0,
       bty = "n")

dev.off()
