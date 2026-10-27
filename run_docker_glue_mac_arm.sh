docker run --rm -it \
  --platform linux/amd64 \
  --name gluetools \
  --network gluetools-network \
  --mount type=bind,source="$PWD",target=/workspace/IFNL-Evolution \
  --workdir /workspace/IFNL-Evolution \
  cvrbioinformatics/gluetools:latest