FROM conda/miniconda3

# docker build -t amcpherson/filtermafs:latest .
# docker push amcpherson/filtermafs:latest

RUN apt-get update
RUN apt-get install git build-essential -y
COPY requirements.txt /requirements.txt
RUN pip install -r requirements.txt

