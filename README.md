# Bayesian-isotope-mixing-model-MixSAIR-for-contribution-rate-calculation
本作品为湖泊面源污染定量溯源完整研究方案，使用 R 语言 MixSIAR 稳定同位素混合模型，耦合碳氮稳定同位素 + 荧光 HIX 三维示踪指标，区分 4 类陆源有机质端元，量化洱海不同入湖河流、旱雨季污染源贡献占比，成果已发表中文核心期刊。提供全流程可复现建模代码、数据输入规范、模型校验绘图、结果输出全套流程，附带季节溯源贡献堆叠柱状结果图，可直接复用至流域有机质、氮磷污染源解析类研究。
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
process_err <- TRUE  # 开启过程误差（水环境溯源标准配置）library(MixSIAR)
getwd()#判断路径，括号里什么都不用填
setwd("C:/Users/Administrator/Documents/2025河流跑代码4端元") #更改路径
mix <- load_mix_data(filename="consumer.csv", 
                     iso_names=c("d13C","d15N","Ratio"), 
                     factors="Sample", 
                     fac_random=FALSE,
                     fac_nested=FALSE,
                     cont_effects=NULL)
mix

source <- load_source_data(filename="source.csv",
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)
source

discr <- load_discr_data(filename="TEF.csv", mix)
discr

plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

plot_prior(alpha.prior=1,source)


model_filename <- "MixSIAR_model.txt"
resid_err <- FALSE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)
jags.1 <- run_model(run="normal", mix, source, discr, model_filename)#步长设置：test/normal/short/long/very long/extreme...

output_JAGS(jags.1, mix, source)

write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# 8. 运行贝叶斯混合模型（normal标准迭代长度，平衡精度与运算速度）
jags.1 <- run_model(
  run = "normal", 
  mix, source, discr, model_filename
)

# 9. 输出模型结果：各采样点、各污染源贡献占比统计+可视化图表
output_JAGS(jags.1, mix, source)
## 运行前置
1. 安装JAGS程序 + R包 `MixSIAR`、`rjags`
2. 修改代码内 `setwd()` 本地文件路径，分步执行即可

## 模型核心设置
- 固定效应分组：采样点Sample
- 开启过程误差、关闭残差误差
- MCMC迭代可选 test/normal/long 精度档位

## 结果规律
1. 旱季农田土壤贡献普遍偏高；雨季冲刷加剧农业源输入，村镇河段生活污水占比上升
2. 山林支流以森林土壤为主，农业河道农田土壤为首要污染源

## 复用说明
替换csv数据即可迁移至湖泊、流域沉积物/氮磷溯源研究
<img width="441" height="185" alt="image" src="https://github.com/user-attachments/assets/262a1ecb-43ef-4409-bb50-49e661e8f6a8" />
<img width="403" height="369" alt="image" src="https://github.com/user-attachments/assets/d6ef7885-0725-48f9-b82e-2f38485ad039" />
