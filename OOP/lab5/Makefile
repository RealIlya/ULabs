TARGET_NAME = lab5

PY_ACTIVATE_BIN = .venv/bin/activate

all: build

run: build
	./$(TARGET_NAME)

build: 
	g++ test/testLab5.cpp ./CauchyEstimate.cpp EmpiricDist.cpp MainDist.cpp MixtureDist.cpp -o $(TARGET_NAME)

py_req:
	chmod +x ./$(PY_ACTIVATE_BIN)
	./$(PY_ACTIVATE_BIN)
	pip install -r requirements.txt

py_plots: py_req build
	python3 ./plot_graphs.py