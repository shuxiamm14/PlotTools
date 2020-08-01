source env.sh

for word in `root-config --cflags`
do
	if [[ "$word" =~ "-std=" ]] ; then
		CXX_STANDARD=`echo $word | tr -dc '0-9'`
		break;
	fi
done

echo "-- Detected root CXX standard:" $CXX_STANDARD

mkdir -p build ; cd build ; cmake .. -DCMAKE_CXX_STANDARD=$CXX_STANDARD ; cd ..