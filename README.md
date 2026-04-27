cal_HPB_auto_iterate.pl 的优势，功能，用法
优势和应用场景：
需要分析的结构较多时，用这个方法，只需要运行一次命令，可避免重复点击鼠标，省时省力
功能：
遍历pdb文件，对每个pdb进行 CDR Annotation（默认IMGT），HPBsidechainSASA calculation，aggregation score calculate，并将相应的结果输出为bsml，csv，dsv。
bsml保存已标注CDR(IMGT)的序列。
csv汇总了所有pdb文件的疏水值，包括HPB，Y，W，HPB+Y+W，相同的数据可在dsv的hierarchy-group看到。
dsv保存了Aggregation propensity Surface,可另存为图片。DSscript 不支持保存图片，只能在脚本运行结束后，逐个打开dsv，另存疏水图。
对于非Fv结构，由于CDR Annotation失效，只能输出整体分子的疏水值。
用法：
第一次使用这脚本，需首先添加环境变量：打开设置界面：右键“此电脑”，选择“属性”，点击“高级系统设置”，然后点击“环境变量”。为Path添加C:\Program Files\BIOVIA\Discovery Studio 2019\bin 
下载pdb结构(通常为ESM batch mode预测的Fv)到本地，层层解压，复制pdb所在文件夹路径(必须纯英文，不能含中文字符)，如'C:\Users\admin\Downloads\struct'。
下载脚本到本地任意自己喜欢的路径。
编辑脚本，修改原$work_dir 为：上一步pdb所在文件夹路径，保存，退出。
在脚本所在文件夹，shift+鼠标右键，打开powershell ，键入 "perl.bat cal_HPB_auto_iterate.pl",回车。等待计算结束，即可在pdb所在文件夹看到dsv bsml csv。
遍历5个结构，大约需要2分钟左右。
