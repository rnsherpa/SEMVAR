FROM python:3.9-slim

WORKDIR /SEMVAR

COPY . /SEMVAR

# Install any necessary dependencies
RUN apt-get update && apt-get install -y \
    && pip install --no-cache-dir -r requirements.txt

# Run SEMVAR.py when the container launches
ENTRYPOINT ["python", "SEMVAR.py"]