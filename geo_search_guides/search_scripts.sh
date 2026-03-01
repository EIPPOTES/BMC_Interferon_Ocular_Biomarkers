
# Bash下载GEO数据的示例
# 下载Series Matrix File
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29801/matrix/GSE29801_series_matrix.txt.gz

# 解压
gunzip GSE29801_series_matrix.txt.gz

# 查看前几行
head -n 20 GSE29801_series_matrix.txt
