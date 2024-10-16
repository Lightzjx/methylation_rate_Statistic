#这个是之间统计DH1c.TE.m.rate的甲基化率的文件内容如下
#chr_1-1024352-1024414   0       3       0.06
#chr_1-1024906-1024948   0       1       0.12
#chr_1-1030940-1031001   0       0       0.0
#chr_1-1037910-1037984   0       2       0.04
#chr_1-1042562-1042638   0       0       0.0
#chr_1-1075531-1075574   0       0       0.0
#脚本原理 分别对2，3列加和，然后相除
import sys
f1 = open(sys.argv[1],"r")
CG,CGM = 0,0
f1.readline()
for i in f1:
	il = i.strip().split()
	CG += int(il[4])
	CGM += int(il[5])
print(f"mean value:\t{CGM/CG}")
