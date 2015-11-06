F=$1
if [ "$1" = "" -o ! -f "$1" ]; then
    echo "$1 doesn't exist."
    exit 1
fi

sed -i "" -e "s/,\ /,/g" ${F}
sed -i "" -e "s/0*,/,/g" ${F}
sed -i "" -e "s/0*\ / /g" ${F}
sed -i "" -e "s/0*)/)/g" ${F}
sed -i "" -e "s/\ ((/((/" ${F}
sed -i "" -e "s/\ (/(/" ${F}
