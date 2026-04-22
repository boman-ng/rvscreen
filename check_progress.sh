while true; do
  size=$(ls -l .local/reference-data/bundle-inputs/host.fa 2>/dev/null | awk '{print $5}')
  echo "Size: $size"
  if ls -l .local/reference-data/bundle-inputs/virus.fa >/dev/null 2>&1; then
    echo "virus.fa appeared!"
    break
  fi
  sleep 60
done
ls -l .local/reference-data/bundle-inputs/
