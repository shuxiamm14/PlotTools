source env.sh

for word in `root-config --cflags`
do
	if [[ "$word" =~ "-std=" ]] ; then
		CXX_STANDARD=`echo $word | tr -dc '0-9'`
		if ((CXX_STANDARD<5)) ; then
			if ((CXX_STANDARD==0)) ; then
				CXX_STANDARD=11;
			elif [[ "$word" =~ "1y" ]] ; then
				CXX_STANDARD=14;
			elif [[ "$word" =~ "1z" ]]; then
				CXX_STANDARD=14;
			else
				echo "ERROR: CXX_STANDARD not found"
				return
			fi
		fi
		break;
	fi
done


echo "-- Detected root CXX standard:" $CXX_STANDARD

mkdir -p build ; cd build ; cmake .. -DCMAKE_CXX_STANDARD=$CXX_STANDARD ; cd ..