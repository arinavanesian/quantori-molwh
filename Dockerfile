FROM python:3.11-slim-bookworm

WORKDIR /app

# Copy requirements from the project `app/` directory (this repository uses `app/` not `src/`)
COPY ./app/requirements.txt .

# Install build deps only for pip compile/build then remove them to keep the image small
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ca-certificates \
    && pip install --no-cache-dir -r requirements.txt \
    && apt-get remove -y --purge build-essential cmake \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /root/.cache/pip

# Copy application package
COPY . .

# Ensure Python can import the `app` package
ENV PYTHONPATH=/app

# Run the FastAPI/uvicorn app from `app.main:app`

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]

# CMD ["python", "-m", "app.main"]