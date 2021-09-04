FROM julia:1.6

WORKDIR /app/
COPY ./src/ /app/

RUN apt-get update
RUN apt-get install -y unzip  # Required for Blink.jl

RUN julia -e 'import Pkg; Pkg.add("Dash")'
RUN julia -e 'import Pkg; Pkg.add("DashCoreComponents")'
RUN julia -e 'import Pkg; Pkg.add("DashHtmlComponents")'
RUN julia -e 'import Pkg; Pkg.add("DataStructures")'
RUN julia -e 'import Pkg; Pkg.add("Distributions")'
RUN julia -e 'import Pkg; Pkg.add("FileIO")'
RUN julia -e 'import Pkg; Pkg.add("ImageMagick")'
RUN julia -e 'import Pkg; Pkg.add("Images")'
RUN julia -e 'import Pkg; Pkg.add("NamedDims")'
RUN julia -e 'import Pkg; Pkg.add("PlotlyJS")'
RUN julia -e 'import Pkg; Pkg.add("YAML")'

EXPOSE 8050
