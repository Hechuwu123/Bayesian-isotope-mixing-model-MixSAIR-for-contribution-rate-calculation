# Bayesian-isotope-mixing-model-MixSAIR-for-contribution-rate-calculation
洱海流域面源污染管控关键问题：区分C3 植物、生活污水、农田土壤、森林土壤4 类有机质端元，量化不同季节、不同入湖采样点污染源贡献，识别旱雨季污染来源差异，为流域控污分区施策提供数据支撑。
# 加载MixSIAR稳定同位素溯源包
library(MixSIAR)
# 1. 工作路径设置（数据文件consumer.csv/source.csv/TEF.csv存放目录）
getwd() # 查看当前工作路径
setwd("C:/Users/Administrator/Documents/") # 切换至数据存储文件夹

# 2. 加载混合样本数据（洱海各采样点有机质样本：δ13C、δ15N、HIX荧光指标多示踪剂）
mix <- load_mix_data(
  filename = "consumer.csv", 
  iso_names = c("d13C","d15N","HIX"), # 碳同位素、氮同位素、荧光HIX三维示踪指标
  factors = "Sample", # 分组变量：各河流采样点
  fac_random = FALSE, # 采样点为固定效应
  fac_nested = FALSE,
  cont_effects = NULL # 无连续环境协变量
)
mix # 查看混合样本数据概览

# 3. 加载4类污染源端元数据（C3植物/生活污水/农田土壤/森林土壤）
source <- load_source_data(
  filename = "source.csv",
  source_factors = NULL, 
  conc_dep = FALSE, # 不考虑浓度依赖校正
  data_type = "means", # 端元输入格式：均值±标准差
  mix = mix
)
source # 查看污染源端元特征数据

# 4. 加载分馏系数TEF数据（同位素营养富集校正值）
discr <- load_discr_data(filename = "TEF.csv", mix = mix)
discr

# 5. 同位素空间分布图绘制（校验端元与样本分离度，溯源有效性检验）
plot_data(
  filename = "isospace_plot", 
  plot_save_pdf = TRUE, 
  plot_save_png = FALSE, 
  mix,source,discr
)

# 6. 先验分布可视化（狄利克雷先验α=1，无偏先验设定）
plot_prior(alpha.prior = 1, source)

# 7. 构建JAGS贝叶斯模型文件
model_filename <- "MixSIAR_model.txt"
resid_err <- FALSE   # 关闭残差误差项
process_err <- TRUE  # 开启过程误差（水环境溯源标准配置）
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# 8. 运行贝叶斯混合模型（normal标准迭代长度，平衡精度与运算速度）
jags.1 <- run_model(
  run = "normal", 
  mix, source, discr, model_filename
)

# 9. 输出模型结果：各采样点、各污染源贡献占比统计+可视化图表
output_JAGS(jags.1, mix, source)
