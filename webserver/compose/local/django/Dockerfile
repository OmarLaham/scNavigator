FROM python:3.8-slim-buster


ENV PYTHONUNBUFFERED 1
ENV PYTHONDONTWRITEBYTECODE 1


ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

# install pkill for celery autoreload
RUN apt-get update && apt-get install -y procps

RUN apt-get update \
  # dependencies for building Python packages
  && apt-get install -y build-essential \
  # psycopg2 dependencies
  && apt-get install -y libpq-dev \
  # Translations dependencies
  && apt-get install -y gettext \
  # cleaning up unused files
  && apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false \
  && rm -rf /var/lib/apt/lists/*

# Requirements are installed here to ensure they will be cached.
COPY ./requirements /requirements

#add conda
RUN apt-get update && apt-get install -y wget \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh
RUN conda --version

#install packages for base conda env
RUN pip install conda-merge==0.2.0 &&\
	conda-merge /requirements/base.yml /requirements/local.yml > /requirements/merged.yml &&\
	pip uninstall --yes conda-merge

RUN conda env update -n root -f /requirements/merged.yml &&\
    conda clean --all --yes

#create separate R single cell env to avoid conflicts
RUN conda env create -f /requirements/r-single-cell-env.yml &&\
    conda clean --all --yes

COPY ./compose/production/django/entrypoint /entrypoint
RUN sed -i 's/\r$//g' /entrypoint
RUN chmod +x /entrypoint

COPY ./compose/local/django/start /start
RUN sed -i 's/\r$//g' /start
RUN chmod +x /start

WORKDIR /app

ENTRYPOINT ["/entrypoint"]
