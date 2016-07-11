FROM ubuntu:14.04
MAINTAINER Steven N Hart, PhD

RUN apt-get update
RUN apt-get install --yes \
 build-essential \
 gcc-multilib \
 apt-utils \
 perl \
 git

# Install perl modules 
RUN apt-get install -y cpanminus

RUN cpanm String::Approx \
	Text::LevenshteinXS

RUN git clone https://github.com/Steven-N-Hart/softsearch.git
RUN cd softsearch && perl install.pl -p $PWD
ENV PERL5LIB=$PERL5LIB:/softsearch/library/