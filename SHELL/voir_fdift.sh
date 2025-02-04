b2a n1=1 <fdift | awk '{print NR " " $1}' > fdift.asc
python ~/PYTHON/plot_xy.py fdift.asc 0 571879 -6. 6. 500



