FROM python:3.9-slim

WORKDIR /SEMVAR

COPY . /SEMVAR

RUN apt-get update && apt-get install -y gosu\
    && pip install --no-cache-dir -r requirements.txt

RUN chmod +x docker-entrypoint.sh

ENTRYPOINT ["./docker-entrypoint.sh"]