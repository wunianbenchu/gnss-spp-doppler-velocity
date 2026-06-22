# gnss-spp-doppler-velocity
C language implementation of GNSS GPS single point positioning (SPP) and Doppler velocity measurement, including RINEX observation file and broadcast ephemeris parsing, least squares positioning solution.
# GNSS SPP 单点定位与多普勒测速
## 项目简介
本项目基于C语言实现了GPS单系统标准单点定位（SPP）与多普勒频移测速算法，完整支持RINEX格式观测数据与广播星历的解析，通过加权最小二乘迭代解算接收机的三维位置与三维速度，可用于卫星导航原理教学验证、实验数据处理与算法二次开发。

## 功能特性
1.  **数据解析**：支持RINEX格式观测文件（.o）与广播星历文件（.p）解析，提取伪距、多普勒观测值及卫星星历参数
2.  **卫星计算**：基于广播星历参数，计算任意历元时刻的GPS卫星三维位置、钟差与钟速
3.  **单点定位**：采用加权最小二乘迭代算法，加入电离层、对流层修正，解算接收机WGS84坐标
4.  **多普勒测速**：利用载波多普勒频移观测值，结合卫星运行速度，解算接收机三维运动速度
5.  **结果输出**：逐历元输出可见星数、位置坐标、速度分量、精度因子（PDOP/HDOP/VDOP）等信息

## 环境依赖
- 编译工具：Visual Studio 2022（推荐，已提供配套工程文件），或任意支持C99标准的C编译器（GCC、Clang等）
- 运行平台：Windows / Linux / macOS
- 输入数据：标准RINEX格式GPS观测文件与广播星历文件

## 文件说明
| 文件名 | 功能说明 |
|--------|----------|
| `main.c` | 核心源码文件，包含数据解析、卫星位置计算、定位测速解算全部逻辑 |
| `brdc1590.24p` | 2024年年积日159的GPS全球广播星历文件 |
| `jfng1590.24o` | jfng测站2024年年积日159的RINEX观测文件 |
| `output.txt` | 程序运行输出的解算结果示例文件 |
| `GNSS单点定位和测速.vcxproj` | Visual Studio 2022 项目配置文件 |
| `GNSS单点定位和测速.vcxproj.filters` | VS项目文件筛选器配置 |

> 注：`x64/` 目录为编译生成的临时产物，无需纳入版本管理

## 编译与运行
### 方式一：Visual Studio 编译（推荐）
1. 下载仓库全部文件，双击 `GNSS单点定位和测速.vcxproj` 用Visual Studio打开
2. 在顶部工具栏选择编译配置（`Release` / `Debug`），平台选择 `x64`
3. 点击「本地Windows调试器」或按快捷键 `F5` 即可编译并运行程序
4. 运行结束后，同目录下的 `output.txt` 即为完整解算结果

### 方式二：GCC 命令行编译
适用于Linux/macOS或Windows MinGW环境：
```bash
# 编译源码，链接数学库
gcc main.c -o gnss_spp_vel -lm

# 运行程序
./gnss_spp_vel
```

## 结果说明
程序输出结果保存在 `output.txt` 中，逐历元包含以下核心信息：
- 观测历元时间（GPS年、月、日、时、分、秒）
- 当前历元可见GPS卫星数量
- 接收机WGS84坐标系三维坐标（X、Y、Z，单位：米）
- 接收机三维速度分量（Vx、Vy、Vz，单位：米/秒）
- 位置精度因子PDOP、平面精度因子HDOP、高程精度因子VDOP

## 注意事项
1. 观测文件与星历文件需与可执行文件置于同一目录；更换测试数据时，请修改 `main.c` 头部的文件名宏定义
2. 本项目默认处理GPS单系统数据，扩展多系统（BDS、GLONASS）可在数据解析模块添加对应逻辑
3. 电离层、对流层采用工程常用简化修正模型，高精度应用可替换为更精密的修正算法

## 参考标准
- GPS 接口控制文件 ICD-GPS-200C
- RINEX 3.x 观测数据与星历数据格式标准
- 卫星导航定位原理与方法
