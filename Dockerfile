# Entorno de la corrida de referencia: Python 3.14.3, semilla de hash fija y un hilo.
FROM python:3.14.3-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=0 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir -e ".[dev]"

COPY . .

CMD ["python", "-m", "pytest", "-q"]
