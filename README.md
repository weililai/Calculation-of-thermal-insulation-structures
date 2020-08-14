# Calculation-of-thermal-insulation-structures
## This is a thermal insulation engineering calculation toolkit
tools based on international standards (ISO 12241 Thermal insulation for building equipment and industrial installations-Calculation rules) are being produced.
## These tools are now included:
### ①Under the requirement of heat loss per unit area or surface temperature, the thickness of pipeline insulation layer is calculated：thickness_of_thermal_Insulation_for_pipe_q_ISO_12241_2008.py

## 这是一个保温工程计算工具包
主要依据 GB50264-2013《工业设备及管道绝热工程设计规范》编写，且有参考 Incropera 的第7版书作《Fundamentals of Heat and Mass Transfer》，也正在开发基于国际标准 ISO 12241 Thermal insulation for building equipment and industrial installations-Calculation rules 的工具包。
## 现在主要包含这些工具：
### ①设定管道保温材料及保温厚度，保温性能计算（这里保温性能主要包括：表面温度、界面温度、热损），并支持多层不同材料复合保温计算：delta_to_temperature_pipe_1.py、delta_to_temperature_pipe_2.py、delta_to_temperature_pipe_3.py。（脚本文件名中的数字表示参与复合的材料种数）
### ②设定管道保温材料、保温厚度及管道长度，保温性能计算（这里保温性能主要包括：表面温度、界面温度、热损、蒸汽温降），并支持多层不同材料复合保温计算：delta_to_temperature_pipe_1_T_change-steam.py、delta_to_temperature_pipe_2_T_change-steam.py、delta_to_temperature_pipe_3_T_change-steam.py。（脚本文件名中的数字表示参与复合的材料种数，该项功能需已安装python库iapws用于解决水蒸气物性计算问题）
### ③给定保温外表面温度及界面温度要求下，保温层厚度计算，并支持多层不同材料复合保温计算：temperature_of_an_compound_insulated_pipe_1.py、temperature_of_an_compound_insulated_pipe_2.py、temperature_of_an_compound_insulated_pipe_3.py。（脚本文件名中的数字表示参与复合的材料种数）
### ④类似的平壁面保温计算工具（文件名中带“wall”的脚本）
## 参考文献：
GB50264-2013《工业设备及管道绝热工程设计规范》
Fundamentals of Heat and Mass Transfer，Incropera 著，第7版
ISO 12241 Thermal insulation for building equipment and industrial installations-Calculation rules
The International Association for the Properties of Water and Steam (IAPWS-IF97)
长输蒸汽管道的温降和压降的计算方法研究_薛永明
蒸汽输热管道的温降设计_赵光显
