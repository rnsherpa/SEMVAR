FROM python:3.9-slim

WORKDIR /SEMVAR

COPY . /SEMVAR

RUN apt-get update && apt-get install -y \
    && pip install --no-cache-dir -r requirements.txt

ENTRYPOINT ["python", "semvar.py"]