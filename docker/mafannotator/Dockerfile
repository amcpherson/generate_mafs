FROM conda/miniconda3

# docker build -t amcpherson/oncokb-annotator:latest .
# docker push amcpherson/oncokb-annotator:latest

RUN apt-get update
RUN apt-get install git -y
RUN git clone --branch raiseonincorrectapikey https://github.com/amcpherson/oncokb-annotator.git
WORKDIR oncokb-annotator
RUN pip install -r requirements/common.txt -r requirements/pip3.txt


