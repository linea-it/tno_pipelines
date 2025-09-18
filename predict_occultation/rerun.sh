
echo "Cleaning up previous run..."
rm -rf /data/outputs/process001

echo "Creating directory and copying files..."
mkdir -p /data/outputs/process001
cp -r /app/predict_occultation/data_sample/asteroid.json /data/outputs/process001
cp -r /app/predict_occultation/data_sample/config.yaml /data/outputs/process001

# cp -r /data/asteroids/2008RH167/2008RH167.bsp /data/outputs/process001
# cp -r /data/asteroids/2008RH167/apmag_and_uncertainties.json /data/outputs/process001

echo "Running prediction..."
/app/predict_occultation/run.sh /data/outputs/process001/config.yaml