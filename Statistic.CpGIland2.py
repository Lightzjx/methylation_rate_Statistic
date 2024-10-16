import sys
f1 = open(sys.argv[1],"r")#~/Temperature/Transcriptome/ref/scaffold/mel6.cpgplot
f2 = open(sys.argv[2],"r")#~/Temperature/methylome/01.Statistic/01.frequency/DH3a.frequency.tsv
islandInfo = {}
for i in f1:
	il = i.strip().split()
	scf = il[0]
	if scf not in islandInfo.keys():
		islandInfo[scf]=[]
	leftplace = il[1]
	rightplace = il[2]
	islandInfo[scf].append([int(leftplace),int(rightplace)])
print("step1")
print(islandInfo.keys())

f2.readline()#去掉第一行
islanddick = {}
for i in f2:
	il = i.split()
	mathyleft,mathyright=int(il[1]),int(il[2])
	try:
		scfisland=islandInfo[il[0]]#[[l,r],[l,r],..,[l,r]]
	except:
		continue
	for islandleft,islandright in scfisland:#遍历上行列表
		islanddickkey = f"{il[0]}-{str(islandleft)}-{str(islandright)}"
		if islanddickkey not in islanddick.keys():
			islanddick[islanddickkey]=[]
		if mathyleft >= islandleft and mathyleft <= islandright:#判断这一小段是否位于island里面，是的话就把这一个行list append加进去
			islanddick[islanddickkey].append(il)
print("step2")
print(len(list(islanddick.keys())))
#print(islanddick)#最终获得全部的island的里面的小碎片的信息dictionary
islandlist = list(islanddick.keys())
islandlist.sort()
islandmathyRate={}
'''for i in islandlist:
	islandmathyInfo = islanddick[i]
	allcg,mathycg = 0,0
	for mathyinfo in islandmathyInfo:
		num_motifs_in_group = int(mathyinfo[0])#小段中的CpG的个数
		methylated_frequency = float(mathyinfo[1])#小段中的CpG甲基化比率
		allcg += num_motifs_in_group#在一个island中的CpG个数持续累加
		if methylated_frequency >= 0.25:#如果甲基化比率>25%,则把这个片段认为是被甲基化的片段，把其中的所有CpG都认为是被甲基化的
			mathycg += num_motifs_in_group#把被甲基化的CpG累加一下
	#islandmathyRate[i]=[mathycg,allcg]
	if allcg==0:#允许在island中没有发现CpG，但是会输出异常
		print("error: this island don't has found any CpG site",i,islandmathyInfo)
		islandmathyRate[i]=[mathycg,allcg,0.00]
	else:
		islandmathyRate[i]=[mathycg,allcg,round(mathycg/allcg,2)]
'''
### 主要函数：片段内mCpG/allCpG   其中：mCpG = Σ(CpGnum*rate); allCpG = Σ(CpGnum); CpGnums为测序最小片段内的CpG的数目，来自tsv文件

#获得两个数据，1、每个CGI中被测到的CpG的个数，2、被测的CpG的甲基化率
for i in islandlist:
	islandmathyInfo = islanddick[i]
	allcg,testedallCpG_nums,testedallCpGm_nums,mathycg,CGms = 0,0,0,0,0
	for mathyinfo in islandmathyInfo:
		num_motifs_in_group = int(mathyinfo[3])#小段中的CpG的个数
		methylated_frequency = float(mathyinfo[6])#小段中的CpG甲基化比率
		CGm = num_motifs_in_group*methylated_frequency
		CGms += CGm
		testedCpG_nums = int(mathyinfo[4])#隔了很久，重新检查，这一行的结果似乎在后面else的第一行以及被注释掉了，所以这一行以及下一行都是不需要的
		testedCpGm_nums = int(mathyinfo[5])

		allcg += num_motifs_in_group#在一个island中的CpG个数持续累加
		testedallCpG_nums += testedCpG_nums
		testedallCpGm_nums += testedCpGm_nums
	if allcg==0:#允许在island中没有发现CpG，但是会输出异常
		print("error: this island don't has found any CpG site",i,islandmathyInfo)
		islandmathyRate[i]=[0,allcg,0.00]
	else:
		#mCpG_rate = round(testedallCpGm_nums/testedallCpG_nums,2)
		mCpG_rate = round(CGms/allcg,2)
		islandmathyRate[i]=[int(CGms),allcg,mCpG_rate]

fw = open(sys.argv[3],"w")#输出文件
def strlist(l):
	newl=[]
	for i in l:
		newl.append(str(i))
	return newl
#CpGi的名称（scf-st-ed）	CpGms(CGms)	allCpGnums(allcg)	mRate(mCpG_rate)	
for i in islandmathyRate.keys():
	il = [i]
	il.extend(islandmathyRate[i])
	il = strlist(il)
	fw.write("\t".join(il)+"\n")
#print(islandmathyRate)
