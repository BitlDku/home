/**
  [ BIT LAB. OFFICIAL ALGORITHM ]   	

  R-value (Biosw site version)
  Ver. 2011-10-19
   
*/
import java.io.*;
import java.util.Arrays; 
import java.util.ArrayList;
import java.util.StringTokenizer;

class calcRvalue 
{
	// PARAMETERS //////////////////////////////////////////////////////////////////////////////

	String dataSetName ;          // class info should be in first column

	int K ;                       // number of nearest Neighbors (default :7) 
	int seta ;                    // thresh hold : determine if a point is belongs to overlapped area. seta <= K (default :3)

	///////////////////////////////////////////////////////////////////////////////////////////
	int NO_OF_DATA ;
	int NO_OF_FEATURE ;
	int NO_OF_CLASS ;

	double [][] totalData ;         
	int [] classInfo ;              
	int [] noClassInstance ;       
	
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	public void findBasicInfo()
	{
		String s1;
		try {
			BufferedReader in = new BufferedReader(new FileReader(dataSetName));
			int cnt = 0, tmpClass =-1, tmpFeature =-1;  
			while ((s1 = in.readLine()) != null)
			{
				StringTokenizer token = new StringTokenizer(s1, ",");
				
				tmpFeature = token.countTokens(); 
				int getClass = Integer.parseInt(token.nextToken().trim());       // class info
				if (getClass > tmpClass)
					tmpClass = getClass;

				cnt++;
			}
			
			in.close(); 

			// basic info -------------------------------------------
			NO_OF_DATA = cnt;
			NO_OF_CLASS = tmpClass + 1;
			NO_OF_FEATURE = tmpFeature -1;	

			System.out.println("## R-value evaluation ##");
			System.out.println("## K="+K+" , seta="+seta);
			System.out.println("Dataset        : " + dataSetName);
			System.out.println("No of instance : " + NO_OF_DATA);
			System.out.println("No of class    : " + NO_OF_CLASS);
			System.out.println("No of feature  : " + NO_OF_FEATURE);
			System.out.println("" );

			totalData = new double[NO_OF_DATA][NO_OF_FEATURE];
			classInfo  = new int[NO_OF_DATA];
			noClassInstance = new int [NO_OF_CLASS];

		} catch (IOException IOe) {
			System.err.println(IOe); // print error message
			System.exit(1);
		}

	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	public void loadData()
	{
		String s1;
		try {
			BufferedReader in = new BufferedReader(new FileReader(dataSetName));
			int i = 0; 
			while ((s1 = in.readLine()) != null)
			{
				StringTokenizer token = new StringTokenizer(s1, ",");

				classInfo[i] = Integer.parseInt(token.nextToken().trim());         // class info
				for (int j=0; j< NO_OF_FEATURE; j++)
				{
					totalData[i][j] = Double.parseDouble(token.nextToken().trim());
				}
				noClassInstance[classInfo[i]]++;

				i++;
			}

			in.close(); 

		// normalize ////////////////////////
		double [] maxVal = new double [NO_OF_FEATURE];

		for (int x=0; x<NO_OF_DATA; x++ )
			for (int j=0; j<NO_OF_FEATURE; j++ )
				if (totalData[x][j] > maxVal[j])
					maxVal[j] = totalData[x][j];

		for (int x=0; x<NO_OF_DATA; x++ )
		{
			for (int j=0; j<NO_OF_FEATURE; j++ ) 
			{
				totalData[x][j] /= maxVal[j];
			}
		}


		} catch (IOException IOe) {
			System.err.println(IOe); // // print error message
			System.exit(1);
		}

	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// Claculate distance between two points  
	public double calcDistance(double [] p1, double [] p2)
	{
		double dist = 0;
		
		for (int i=0; i<NO_OF_FEATURE; i++)
			dist += Math.pow(p1[i]-p2[i],2);
			//dist = Math.sqrt(dist);   --> can be omitted
		return dist;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// find k-nearest neighbors 
	
	public ArrayList<Integer> findKnn(double [] p1, int p1_idx, double [][] targetList)
	{
		ArrayList<Double> winner = new ArrayList<Double>(K);
		ArrayList<Integer> winnerIdx = new ArrayList<Integer>(K);

		for (int j=0; j<targetList.length;j++)
		{
			if (p1_idx == j) continue ;     // skip if instance (i) meets itself 

			double tmpDistance = calcDistance (p1, targetList[j]);

			int insertPosition = winner.size();
			for (int x=winner.size()-1; x >= 0; x--)
				if (tmpDistance > winner.get(x).doubleValue())
					break ;
				else
					insertPosition = x;
			
			winner.add(insertPosition, tmpDistance);                      // distance
			winnerIdx.add(insertPosition, j);                             // class info

			if (winner.size() > K)
			{   
				winner.remove(K); 
				winnerIdx.remove(K); 
			}
		}

		return winnerIdx ; 
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// R-value : R(Ci,Cj)  
	public double subClacRvalue1(double[][] arr1, double[][] arr2)
	{
		int [] Rvalue1 = new int [arr1.length] ;
		double [][] totalArr = new double[arr1.length+arr2.length][NO_OF_FEATURE+1] ; 	

		// megrge arr1, arr2 -------------------------- 
		for (int i=0;i<arr1.length; i++ )
		{
			for (int k=0; k<NO_OF_FEATURE; k++ )
				totalArr[i][k] = arr1[i][k];

			totalArr[i][NO_OF_FEATURE] = 1;
		}
		for (int i=0;i<arr2.length; i++ )
		{
			for (int k=0; k<NO_OF_FEATURE; k++ )
				totalArr[arr1.length+i][k] = arr2[i][k];

			totalArr[arr1.length+i][NO_OF_FEATURE] = 2;
		}

		// calc KNN -----------------------------------
		for (int i=0;i<arr1.length;i++) 
		{
			ArrayList<Integer> winnerIdx ;

		    winnerIdx = findKnn(arr1[i], i, totalArr);

			for (int x=0; x<winnerIdx.size(); x++)
				if (totalArr[winnerIdx.get(x).intValue()][NO_OF_FEATURE] != 1)  // 1 : class of arr1
					Rvalue1[i]++;
		}


		for (int k=0; k<Rvalue1.length; k++)
			if (Rvalue1[k] > seta)
				Rvalue1[k] = 1;
		    else
				Rvalue1[k] = 0;

		double avg1=0 ;
		for (int j=0;j<arr1.length;j++)
			avg1 += Rvalue1[j];

		return avg1;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// R-value : R(Ci)  
	public double subClacRvalue2(double[][] arr1, int cIndex)
	{
		int K = 7;
		int [] Rvalue1 = new int [arr1.length] ;

		// calc KNN -----------------------------------
		for (int i=0;i<arr1.length;i++) 
		{
			ArrayList<Integer> winnerIdx ;

		    winnerIdx = findKnn(arr1[i], (int)arr1[i][NO_OF_FEATURE], totalData);

			for (int x=0; x<winnerIdx.size(); x++)
				if (classInfo[winnerIdx.get(x).intValue()] != cIndex)   
					Rvalue1[i]++;
		}

		for (int k=0; k<Rvalue1.length; k++)
		{
			if (Rvalue1[k] > seta)
				Rvalue1[k] = 1;
		    else
				Rvalue1[k] = 0;
		}

		double avg1=0 ;
		for (int j=0;j<arr1.length;j++)
		{
			avg1 += Rvalue1[j];
		}

		return avg1;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////	
	public void calcRvalue()
	{

		try {
			String outputData = dataSetName + "_rvalue.txt";
			FileWriter   fw_1 = new FileWriter(outputData);  
			BufferedWriter bw_1 = new BufferedWriter(fw_1);   
			PrintWriter outFile_1 = new PrintWriter(bw_1); 

			outFile_1.println("** " + dataSetName); 
			outFile_1.println();
			// calculate R-value (R(C1,C2)-----------------------------------------
			for (int i=0; i<NO_OF_CLASS; i++)
			{
				double [][] arr1 ; //= new double [100][2];		
				double [][] arr2 ; //= new double [100][2];	
				for (int j=0; j<NO_OF_CLASS; j++)
					if (i < j)
					{
						// get arr1, array2
						arr1 = new double[noClassInstance[i]][NO_OF_FEATURE];
						arr2 = new double[noClassInstance[j]][NO_OF_FEATURE];
						int cnt1=0, cnt2=0;

						for (int k=0; k<NO_OF_DATA; k++)
						{
							// arr1
							if (classInfo[k] == i)
							{
								for (int x=0; x<NO_OF_FEATURE; x++)
									arr1[cnt1][x] = totalData[k][x];
								cnt1++;
							}
							// arr2
							if (classInfo[k] == j)
							{
								for (int x=0; x<NO_OF_FEATURE; x++)
									arr2[cnt2][x] = totalData[k][x];
								cnt2++;
							}
						}
						// call calculation function
						String cGroup = "("+String.valueOf(i+1)+","+String.valueOf(j+1)+")";
						double tmpSum = subClacRvalue1(arr1, arr2) +
										subClacRvalue1(arr2, arr1) ;
						tmpSum = (tmpSum)/(arr1.length + arr2.length);
						System.out.println("R "+ cGroup + " = " + tmpSum);
						outFile_1.println("R "+ cGroup + " = " + tmpSum);
					}
			}

			// calculate R-value (R(C1),R(f))-----------------------------------------
			double RC1 =0, Rf =0;
			for (int i=0; i<NO_OF_CLASS; i++)
			{
				double [][] arr1 = new double [noClassInstance[i]][NO_OF_FEATURE+1];		  // 
				RC1 =0;                                       // initialize
				int cnt1 =0;
				// get arr1
				for (int k=0; k<NO_OF_DATA; k++)
					if (classInfo[k] == i)
					{
						for (int x=0; x<NO_OF_FEATURE; x++)
							arr1[cnt1][x] = totalData[k][x];
						arr1[cnt1][NO_OF_FEATURE] = k;          // location info. 

						cnt1++;
					}

				// call calculation function
				String cGroup = "R("+String.valueOf(i+1)+")";
				RC1 = subClacRvalue2(arr1, i);         // R-value R(C1)

				Rf += RC1; 
				RC1 /= (double)arr1.length ;    // 
				System.out.println(cGroup + " = " + RC1);
				outFile_1.println(cGroup + " = " + RC1);

			}

			// Print R(f)
			Rf /= (double)totalData.length ;    // 
			System.out.println("R(f)= "+ Rf);
			outFile_1.println("R(f)= "+ Rf);

			outFile_1.close();


		} catch (IOException IOe) {
			System.err.println(IOe);
			System.exit(1);
		}	

	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	public static void main(String[] args) 
	{
		calcRvalue obj = new calcRvalue();
		boolean error = false ;
		// argument parsing
		switch (args.length) {
			case 0: System.out.println("(Error) No dataset file.");
					break;
			case 1: obj.dataSetName = args[0]; obj.K = 7; obj.seta = 3 ; 
                    break;
			case 3: obj.dataSetName = args[0]; 
		            obj.K = Integer.parseInt(args[1]); 
				    obj.seta = Integer.parseInt(args[2]) ;
                    break;
			default: System.out.println("Wrong parameters!");
					error = true;
                    break;
		}		

		if (!error)
		{
			obj.findBasicInfo();
			obj.loadData();
			obj.calcRvalue();
		}
		
	}
}
