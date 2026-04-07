FOREFIRE_IMAGE := forefire
APP_PORT       := 8501
REPO_ROOT      := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../..)

.PHONY: app build help

help:
	@echo "make build  — build the forefire Docker image (run once)"
	@echo "make app    — launch the Streamlit app on http://localhost:$(APP_PORT)"

## Build the underlying forefire Docker image (required by the app)
build:
	docker build -t $(FOREFIRE_IMAGE) $(REPO_ROOT)

## Launch the Streamlit app on the host. Requires `docker` and `streamlit`/`folium`
## installed locally; the app shells out to `docker run forefire ...` per simulation.
app:
	@echo "Starting ForeFire Spread Predictor on http://localhost:$(APP_PORT)"
	streamlit run $(dir $(lastword $(MAKEFILE_LIST)))streamlit_app.py --server.port=$(APP_PORT)
