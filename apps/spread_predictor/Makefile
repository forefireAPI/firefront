FOREFIRE_IMAGE := forefire
APP_PORT       := 8501
REPO_ROOT      := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../..)

.PHONY: app build ensure-image help

help:
	@echo "make build  — (re)build the forefire Docker image"
	@echo "make app    — build the image if missing, then launch the Streamlit app"

## (Re)build the underlying forefire Docker image
build:
	docker build -t $(FOREFIRE_IMAGE) $(REPO_ROOT)

## Build the forefire image only if it doesn't already exist
ensure-image:
	@if [ -z "$$(docker images -q $(FOREFIRE_IMAGE) 2>/dev/null)" ]; then \
		echo "forefire image not found — building it now..."; \
		$(MAKE) build; \
	else \
		echo "forefire image already present."; \
	fi

## Launch the Streamlit app on the host. Requires `docker` and `streamlit`/`folium`
## installed locally; the app shells out to `docker run forefire ...` per simulation.
app: ensure-image
	@echo "Starting ForeFire Spread Predictor on http://localhost:$(APP_PORT)"
	streamlit run $(dir $(lastword $(MAKEFILE_LIST)))streamlit_app.py --server.port=$(APP_PORT)
