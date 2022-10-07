
FROM python:3.8-slim-buster

ENV PYTHONUNBUFFERED 1


ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

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

RUN addgroup --system django \
    && adduser --system --ingroup django django



# Requirements are installed here to ensure they will be cached.
COPY ./requirements /requirements

#add conda
RUN apt-get update && apt-get install -y wget \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

RUN pip install conda-merge==0.2.0 &&\
	conda-merge /requirements/base.yml /requirements/local.yml > /requirements/merged.yml &&\
	pip uninstall --yes conda-merge

RUN conda env update -n root -f /requirements/merged.yml &&\
    conda clean --all --yes
    rm -rf /requirements

COPY --chown=django:django ./compose/production/django/entrypoint /entrypoint
RUN sed -i 's/\r$//g' /entrypoint
RUN chmod +x /entrypoint


COPY --chown=django:django ./compose/production/django/start /start
RUN sed -i 's/\r$//g' /start
RUN chmod +x /start


COPY --chown=django ./compose/production/django/create_superuser /create_superuser
RUN sed -i 's/\r//' /create_superuser
RUN chmod +x /create_superuser

RUN mkdir -p /app && chown django /app

COPY --chown=django . /app

RUN mkdir /app/staticfiles && chown -R django /app/staticfiles
#mkdir -p means make dir only if doesn't exist
RUN mkdir -p /app/planet/media && chown -R django /app/planet/media

USER django

WORKDIR /app

ENTRYPOINT ["/entrypoint"]
