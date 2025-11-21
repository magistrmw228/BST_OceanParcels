#!/bin/bash

URL="https://downloader.disk.yandex.ru/disk/72f3aa585c2274c2099e833f43df754280bf36615ab5b4e9fc24c800e8a1233d/6920cc40/fKqInKw3d7bLFOeFnMGnhLBQUiXTWVgy20YsbQLq54VHelyNtGXu05CToott-ljiG6SLj-zxG6q46muDJddL0P91y7_2RP5S3igAtfWXUkqr8npumZHI4midPdWhecNq?uid=1130000069824731&filename=Input.zip&disposition=attachment&hash=&limit=0&content_type=application%2Fzip&owner_uid=1130000069824731&fsize=1466312394&hid=f27ee31ad41830a697465b7a65587d1a&media_type=compressed&tknv=v3&etag=3ec2f6c24b64cf995d93756518cdb494"
ZIP_FILE="Input.zip"

echo "Downloading files (1.4 Gb need to be downloaded)..."
wget -O "$ZIP_FILE" "$URL"

if [ $? -ne 0 ]; then
    echo "Error of downloading"
    exit 1
fi

echo "Unzip files..."
unzip -o "$ZIP_FILE"

if [ $? -ne 0 ]; then
    echo "Unzip error"
    exit 1
fi

echo "Done! Input files are in Input/ folder"
rm "$ZIP_FILE"
