APP_IMAGE := forefire-streamlit
APP_PORT  := 8501

.PHONY: app app-build app-stop stop test

## Build the Streamlit app Docker image
app-build:
	docker build -f Dockerfile.streamlit -t $(APP_IMAGE) .

## Run the Streamlit app (stops previous, builds, runs)
app: app-stop app-build
	@echo "Starting ForeFire Streamlit app on http://localhost:$(APP_PORT)"
	docker run --rm -p $(APP_PORT):8501 --name $(APP_IMAGE) $(APP_IMAGE)

## Stop the running app
app-stop:
	@docker stop $(APP_IMAGE) 2>/dev/null || true
	@docker rm $(APP_IMAGE) 2>/dev/null || true

## Stop the running app (alias)
stop: app-stop

## Run CI tests in Docker
test:
	docker build -t forefire-test .
	docker run --rm forefire-test bash -c "cd tests/runff && bash ff-run.bash"
