/*
 * HashTable.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#include "HashTable.h"

DataVector::DataVector()
{
	dataList = NULL;
	keystring = "";
}
DataVector::~DataVector()
{
	if(dataList!=NULL)
	{
		dataList->clear();
		delete dataList;
	}

}
HashTable::HashTable() {
	// TODO Auto-generated constructor stub
	hashKeyLength = Config::hashKeyLength;
	hashTableName = "";
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
	hashTableSize = 0;
	hashTable = NULL;
}
HashTable::HashTable(UINT16 keylength) {
	// TODO Auto-generated constructor stub
	hashKeyLength = keylength;
	hashTableName = "";
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
	hashTableSize = 0;
	hashTable = NULL;

}
HashTable::HashTable(UINT16 keylength, string name) {
	// TODO Auto-generated constructor stub
	hashKeyLength = keylength;
	hashTableName = name;
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
	hashTableSize = 0;
	hashTable = NULL;

}
HashTable::HashTable(UINT16 keylength, UINT64 dataSetSize) {
	// TODO Auto-generated constructor stub
	hashKeyLength = keylength;
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
	InitializeWithDataSize(dataSetSize);

}

UINT64 HashTable::getMaximumHashTableCollision()
{
	return this->maxSingleHashCollision;
}
void HashTable::InitializeWithDataSize(UINT64 dataSetSize)
{
	UINT64 size = getPrimeLargerThanNumber(dataSetSize * 2 + 1);  // Size should be at least twice the number of entries in the hash table to reduce hash collision.
	setHashTableSizeAndInitialize(size);
}

HashTable::~HashTable() {
	// TODO Auto-generated destructor stub
	if(hashTable!=NULL)
	{
		for(UINT64 i=0; i<hashTable->size();i++)
			delete hashTable->at(i);
		hashTable->clear();
		delete hashTable;
	}
}

//This hash function only takes care of the first 32 bps and the last 32 bps of the total key length.
//What's in between will not be hashed.
//So the in between sequence will be compared after the key is hit in the hash table, and if the in between sequence
//doesn't match this method will report a collision.
UINT64 HashTable::hashFunction(const string & subString)
{
	UINT64 sum1 = 1, sum2 = 1, length = subString.length();// initial value of 1 to avoid zero when doing the multiplication.
	for(UINT64 i = 0; i < length; i++)			// We take the bit representation of the string. A = 00, C = 01, G = 11 and T = 10
	{
		if(i < 32)
			// sum1 is for the first 32 bp. bit shifted to left 2 bits.
			// Change the character to integer. A=Ox41=01000001
			//                                  C=0x43=01000011
			//                                  G=0x47=01000111
			//                                  T=0x54=01010100
			// Then, shift to right way 1 bit.
			// Then, bit and operation with 00000011
			// Then, it just have A=00, C=01, G=11,T=10
			sum1 = (sum1 << 2) | ( ( (int)(subString[i] ) >> 1 ) & 0X03 );
		else
			sum2 = (sum2 << 2) | ( ( (int)(subString[i] ) >> 1 ) & 0X03 );
	}
	// multiply two mod results, and mod again.
	return ((sum1 % hashTableSize) * (sum2  % hashTableSize)) % hashTableSize; 	// Modulus operation to get the index in the hash table.
}

//Function to get a prime number larger than a given number.
//For efficiency of the hash function the hash table size is set to a prime number.
//This will reduce the number of hash collision.
UINT64 HashTable::getPrimeLargerThanNumber(UINT64 number)
{
	// Pre-computed list of prime number.
	UINT64 array[450] = {1114523, 1180043, 1245227, 1310759, 1376447, 1442087, 1507379, 1573667, 1638899, 1704023, 1769627, 1835027, 1900667, 1966127, 2031839, 2228483, 2359559, 2490707, 2621447, 2752679, 2883767, 3015527, 3145739, 3277283, 3408323, 3539267, 3670259, 3801143, 3932483, 4063559, 4456643, 4718699, 4980827, 5243003, 5505239, 5767187, 6029603, 6291563, 6553979, 6816527, 7079159, 7340639, 7602359, 7864799, 8126747, 8913119, 9437399, 9962207, 10485767, 11010383, 11534819, 12059123, 12583007, 13107923, 13631819, 14156543, 14680067, 15204467, 15729647, 16253423, 17825999, 18874379, 19923227, 20971799, 22020227, 23069447, 24117683, 25166423, 26214743, 27264047, 28312007, 29360147, 30410483, 31457627, 32505983, 35651783, 37749983, 39845987, 41943347, 44040383, 46137887, 48234623, 50331707, 52429067, 54526019, 56623367, 58720307, 60817763, 62915459, 65012279, 71303567, 75497999, 79691867, 83886983, 88080527, 92275307, 96470447, 100663439, 104858387, 109052183, 113246699, 117440699, 121635467, 125829239, 130023683, 142606379, 150994979, 159383759, 167772239, 176160779, 184549559, 192938003, 201327359, 209715719, 218104427, 226493747, 234882239, 243269639, 251659139, 260047367, 285215507, 301989959, 318767927, 335544323, 352321643, 369100463, 385876703, 402654059, 419432243, 436208447, 452986103, 469762067, 486539519, 503316623, 520094747, 570425399, 603979919, 637534763, 671089283, 704643287, 738198347, 771752363, 805307963, 838861103, 872415239, 905971007, 939525143, 973079279, 1006633283, 1040187419, 1140852767, 1207960679, 1275069143, 1342177379, 1409288183, 1476395699, 1543504343, 1610613119, 1677721667, 1744830587, 1811940419, 1879049087, 1946157419, 2013265967, 2080375127, 2281701827, 2415920939, 2550137039, 2684355383, 2818572539, 2952791147, 3087008663, 3221226167, 3355444187, 3489661079, 3623878823, 3758096939, 3892314659, 4026532187, 4160749883, 4563403379, 4831838783, 5100273923, 5368709219, 5637144743, 5905580687, 6174015503, 6442452119, 6710886467, 6979322123, 7247758307, 7516193123, 7784629079, 8053065599, 8321499203, 9126806147, 9663676523, 10200548819, 10737418883, 11274289319, 11811160139, 12348031523, 12884902223, 13421772839, 13958645543, 14495515943, 15032386163, 15569257247, 16106127887, 16642998803, 18253612127, 19327353083, 20401094843, 21474837719, 22548578579, 23622320927, 24696062387, 25769803799, 26843546243, 27917287907, 28991030759, 30064772327, 31138513067, 32212254947, 33285996803, 36507222923, 38654706323, 40802189423, 42949673423, 45097157927, 47244640319, 49392124247, 51539607599, 53687092307, 55834576979, 57982058579, 60129542339, 62277026327, 64424509847, 66571993199, 73014444299, 77309412407, 81604379243, 85899346727, 90194314103, 94489281203, 98784255863, 103079215439, 107374183703, 111669150239, 115964117999, 120259085183, 124554051983, 128849019059, 133143986399, 146028888179, 154618823603, 163208757527, 171798693719, 180388628579, 188978561207, 197568495647, 206158430447, 214748365067, 223338303719, 231928234787, 240518168603, 249108103547, 257698038539, 266287975727, 292057776239, 309237645803, 326417515547, 343597385507, 360777253763, 377957124803, 395136991499, 412316861267, 429496730879, 446676599987, 463856468987, 481036337207, 498216206387, 515396078039, 532575944723, 584115552323, 618475290887, 652835029643, 687194768879, 721554506879, 755914244627, 790273985219, 824633721383, 858993459587, 893353198763, 927712936643, 962072674643, 996432414899, 1030792152539, 1065151889507, 1168231105859, 1236950582039, 1305670059983, 1374389535587, 1443109012607, 1511828491883, 1580547965639, 1649267441747, 1717986918839, 1786706397767, 1855425872459, 1924145348627, 1992864827099, 2061584304323, 2130303780503, 2336462210183, 2473901164367, 2611340118887, 2748779070239, 2886218024939, 3023656976507, 3161095931639, 3298534883999, 3435973836983, 3573412791647, 3710851743923, 3848290698467, 3985729653707, 4123168604483, 4260607557707, 4672924419707, 4947802331663, 5222680234139, 5497558138979, 5772436047947, 6047313952943, 6322191860339, 6597069767699, 6871947674003, 7146825580703, 7421703488567, 7696581395627, 7971459304163, 8246337210659, 8521215117407, 9345848837267, 9895604651243, 10445360463947, 10995116279639, 11544872100683, 12094627906847, 12644383722779, 13194139536659, 13743895350023, 14293651161443, 14843406975659, 15393162789503, 15942918604343, 16492674420863, 17042430234443, 18691697672867, 19791209300867, 20890720927823, 21990232555703, 23089744183799, 24189255814847, 25288767440099, 26388279068903, 27487790694887, 28587302323787, 29686813951463, 30786325577867, 31885837205567, 32985348833687, 34084860462083, 37383395344739, 39582418600883, 41781441856823, 43980465111383, 46179488367203, 48378511622303, 50577534878987, 52776558134423, 54975581392583, 57174604644503, 59373627900407, 61572651156383, 63771674412287, 65970697666967, 68169720924167, 74766790688867, 79164837200927, 83562883712027, 87960930223163, 92358976733483, 96757023247427, 101155069756823, 105553116266999, 109951162779203, 114349209290003, 118747255800179, 123145302311783, 127543348823027, 131941395333479, 136339441846019, 149533581378263, 158329674402959, 167125767424739, 175921860444599, 184717953466703, 193514046490343, 202310139514283, 211106232536699, 219902325558107, 228698418578879, 237494511600287, 246290604623279, 255086697645023, 263882790666959, 272678883689987, 299067162755363, 316659348799919, 334251534845303, 351843720890723, 369435906934019, 387028092977819, 404620279022447, 422212465067447, 439804651111103, 457396837157483, 474989023199423, 492581209246163, 510173395291199, 527765581341227, 545357767379483, 598134325510343, 633318697599023, 668503069688723, 703687441776707, 738871813866287, 774056185954967, 809240558043419, 844424930134187, 879609302222207, 914793674313899, 949978046398607, 985162418489267, 1020346790579903, 1055531162666507, 1090715534754863}; 	//list of 450 sorted prime numbers;
	vector<UINT64> primeNumbers(array, array + sizeof array / sizeof array[0]);
	for(UINT64 i = 0; i < primeNumbers.size(); i++)
		if(primeNumbers.at(i)> number)
			return primeNumbers.at(i+1);		// Return the smallest prime in the list larger than number.
	return number + 1;
}

void HashTable::setHashTableSizeAndInitialize(UINT64 size)
{
	cout << "Hash Table size set to: " << size << endl;
	hashTableSize=size;

	hashTable = new vector < DataVector *>();
	hashTable->reserve(size);

	for(UINT64 i = 0; i < hashTable->capacity(); i++) // Initialize the hash table.
	{
		vector<UINT64> * newList = new vector<UINT64>();
		DataVector * dataVector = new DataVector();
		dataVector->dataList = newList;
		newList->resize(newList->size());
		hashTable->push_back(dataVector);
	}
}


bool HashTable::findInsertIndex(string keyString, UINT64& hashTableIndex)
{
	UINT64 currentCollisionNum = 0;

	hashTableIndex = hashFunction(keyString);
	bool isEmpty = hashTable->at(hashTableIndex)->dataList->empty();
	while(!isEmpty)
	{
		if(currentCollisionNum>1000 || currentCollisionNum>hashTableSize)
		{
			cout<<"Too much collision, need to enlarge the hashTableSize"<<endl;
			return false;
		}
		if(keyString == hashTable->at(hashTableIndex)->keystring)
				break;
		numberOfHashCollision++;
		currentCollisionNum++;
		hashTableIndex = (hashTableIndex == hashTableSize - 1) ? 0: hashTableIndex + 1; 	// Increment the index
		isEmpty = hashTable->at(hashTableIndex)->dataList->empty();
	}

	if(currentCollisionNum> this->maxSingleHashCollision)
		this->maxSingleHashCollision = currentCollisionNum;

	return true;
}
bool HashTable::insertIntoHashTable(string keyString, UINT64 readID)
{

	UINT64 hashTableIndex;
	bool insertSuccess = findInsertIndex(keyString, hashTableIndex);


	hashTable->at(hashTableIndex)->dataList->push_back(readID);							// Add the string in the list.
	hashTable->at(hashTableIndex)->dataList->resize(hashTable->at(hashTableIndex)->dataList->size());		// Resize to reduce space.
	if(hashTable->at(hashTableIndex)->keystring=="")
		hashTable->at(hashTableIndex)->keystring = keyString;
	else if(hashTable->at(hashTableIndex)->keystring!=keyString)
	{
		cout<<"contradiction of inserting keystring in function insertIntoHashTable()"<<endl;
		return false;
	}

	return insertSuccess;

}

bool HashTable::getIndexFromHashTable(string subString, UINT64& hashTableIndex)
{
	UINT64 currentCollisionNum = 0;
	hashTableIndex = hashFunction(subString);
	vector<UINT64> tempvector = (*hashTable->at(hashTableIndex)->dataList);
	bool isEmpty = hashTable->at(hashTableIndex)->dataList->empty();
	while(!isEmpty && currentCollisionNum<=maxSingleHashCollision)
	{

			string keysequence = hashTable->at(hashTableIndex)->keystring;
			if(subString == keysequence)
					return true;

			currentCollisionNum++;
			hashTableIndex = (hashTableIndex == hashTableSize - 1) ? 0: hashTableIndex + 1; 	// Increment the index
			isEmpty = hashTable->at(hashTableIndex)->dataList->empty();
	}
	return false;
}
string HashTable::getHashKeyFromIndex(UINT64 hashTableIndex)
{
	string keyString = hashTable->at(hashTableIndex)->keystring;
	return keyString;
}
vector<UINT64> * HashTable::getReadIDListOfReads(string subString)
{
	UINT64 hashTableIndex;
	if(getIndexFromHashTable(subString,hashTableIndex))
		return hashTable->at(hashTableIndex)->dataList;
	else return NULL;
}