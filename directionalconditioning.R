
df <- read.table(text = "
phenotype                                Tukey_Z2       pval     n_shared    n_total    n_phenotype
'Age at which the menopause started'    -0.14595278    1.000    476         18121       18681
'Chronic insomnia'                       0.71341236    0.023    14          198         202
'Daily cigarrete'                         -0.52246419    1.000    72          2862        3062
'Coronary artery disease'               -0.70206204    1.000    306         11717       12173
'Covid-19 infection'                     0.22538920    0.086    35          1323        1390
'Covid-19 Very severe'                 -0.45984234    0.999    59          3052        3138
'Depressive disorder'                  -0.05306575    0.761    88          3384        3566
'Diabetes mellitus type 1'             -0.45067544    1.000    144         17445       18974
'Diabetes mellitus type 2'              0.09272861    0.032    314         10052       12359
'Diet'                                  0.57301349    0.000    207         9148        10987
'Drink Weekly'                         -0.10198629    0.948    164         6533        7546
'Female Infertility'                   -1.42738085    0.991    7           61          68
'Heart failure'                        -1.31665554    1.000    31          513         538
'Hypertension'                         -0.52365223    1.000    470         14842       16278
'Influenza infection'                  -0.63108917    0.862    2           60          60
'Lewy body dementia'                   71.21064779    0.000    5           68          128
'Loneliness'                            0.20498915    0.247    9           176         181
'Obesity'                              -2.51293603    1.000    6           329         337
'Parkinson disease'                     2.76595919    0.000    45          2534        2753
'Psychosocial stress measurement'      -0.60742467    0.786    1           2           2
'Skin color'                           -0.66084933    0.800    3           3           14
'Sleep apnea'                          -0.65178939    1.000    39          772         793
'Stroke'                               -1.95512113    1.000    19          367         375
", header = TRUE)


library(ggplot2)

# Order phenotypes by descending Z² (Lewy body dementia on top)
df$phenotype <- factor(df$phenotype, levels = df$phenotype[order(-df$Tukey_Z2)])

# Define significance based on abs(1 - pval) < 0.05
df$significant <- abs(df$pval - 1) > 0.05 

# Plot
ggplot(df, aes(x = phenotype, y = Tukey_Z2, fill = significant)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "gray80"),
                    name = "p-value",
                    labels = c("Not significant", "Significant")) +
  labs(
    title = "Directional Conditioning Between Phenotypes",
    x = "",
    y = "Tukey's Mean Z-score²"
  ) +
  theme_minimal(base_size = 12)
