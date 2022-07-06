
# 1st Psub/Pwell scan Goe DUT
#python3 tj2-reco.py --rawfile /home/benjamin/desy_tb/run279.txt --runno 279
#python3 tj2-reco.py --rawfile /home/benjamin/desy_tb/run280.txt --runno 280
#python3 tj2-reco.py --rawfile /home/benjamin/desy_tb/run281.txt --runno 281
#python3 tj2-reco.py --rawfile /home/benjamin/desy_tb/run282.txt --runno 282
#python3 tj2-reco.py --rawfile /home/benjamin/desy_tb/run283.txt --runno 283

python3 histo-plotter-tj2.py --ifile=root-files/Histos-TJ2-run279-reco.root
python3 histo-plotter-tj2.py --ifile=root-files/Histos-TJ2-run280-reco.root
python3 histo-plotter-tj2.py --ifile=root-files/Histos-TJ2-run281-reco.root
python3 histo-plotter-tj2.py --ifile=root-files/Histos-TJ2-run282-reco.root
python3 histo-plotter-tj2.py --ifile=root-files/Histos-TJ2-run283-reco.root

# 2nd Psub/Pwell scan Goe DUT
python3 tj2-reco.py --rawfile /home/benjamin/desy_tb/run377.txt --runno 377
python3 tj2-reco.py --rawfile /home/benjamin/desy_tb/run378.txt --runno 378
python3 tj2-reco.py --rawfile /home/benjamin/desy_tb/run379.txt --runno 379
python3 tj2-reco.py --rawfile /home/benjamin/desy_tb/run380.txt --runno 380

python3 histo-plotter-tj2.py --ifile=root-files/Histos-TJ2-run377-reco.root
python3 histo-plotter-tj2.py --ifile=root-files/Histos-TJ2-run378-reco.root
python3 histo-plotter-tj2.py --ifile=root-files/Histos-TJ2-run379-reco.root
python3 histo-plotter-tj2.py --ifile=root-files/Histos-TJ2-run380-reco.root



