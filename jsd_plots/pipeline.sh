cp rates.csv rates.csv.bak

# throw away iron oxidation on oxygen
drop_column.pl 0 rates.csv.bak > rates.csv

# convert the rates rows into a matrix of JSDs
python -c "import survey.analyze as sa; sa.csv2jsd()"
