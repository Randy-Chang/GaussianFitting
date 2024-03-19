/*
 * 參考文獻
 * COMPARISON OF ALGORITHMS FOR FITTING A GAUSSIAN FUNCTION USED IN TESTING SMART SENSORS
 * 網址 https://intapi.sciendo.com/pdf/10.2478/jee-2015-0029
*/

namespace GaussianFitting
{
    enum EMethod { TwoPoint, CARUANAS_ALGORITHM, GUOS_ALGORITHM }
    public partial class Form1 : Form
    {
        FarFieldData FFData;
        GaussianParameter GParam;
        List<double> yFFVG, yFFHG;
        EMethod method;

        public Form1()
        {
            InitializeComponent();

            string dataPath = @"C:\Users\gdba5\Desktop\Input.csv";
            string dataPathSave = @"C:\Users\gdba5\Desktop\Output.csv";

            LoadData(dataPath, out FFData);

            method = EMethod.GUOS_ALGORITHM;

            switch (method)
            {
                case EMethod.TwoPoint:
                    InputConditionUsingTwoPoint(FFData.xFFV, FFData.yFFV, 50, out GParam);
                    CalculateGaussian(FFData.xFFV, GParam, out yFFVG);
                    break;

                case EMethod.CARUANAS_ALGORITHM:
                    CARUANAS_ALGORITHM(FFData.xFFV, FFData.yFFV, out GParam);
                    CalculateGaussian(FFData.xFFV, GParam, out yFFVG);
                    break;

                case EMethod.GUOS_ALGORITHM:
                    GUOS_ALGORITHM(FFData.xFFV, FFData.yFFV, out GParam);
                    CalculateGaussian(FFData.xFFV, GParam, out yFFVG);
                    break;

                default: break;
            }
           

            SaveData(dataPathSave, FFData, yFFHG, yFFVG);
        }

        #region Gaussian Fittin CARUANAS_ALGORITHM

        void CARUANAS_ALGORITHM(List<double> xData, List<double> yData,
                                    out GaussianParameter gaussianParameter)
        {
            gaussianParameter = new GaussianParameter();

            double term11 = 0;
            double term12 = 0;
            double term13 = 0;
            double term21 = 0;
            double term22 = 0;
            double term23 = 0;
            double term31 = 0;
            double term32 = 0;
            double term33 = 0;
            double termY1 = 0;
            double termY2 = 0;
            double termY3 = 0;

            int count = yData.Count;

            for (int i = 0; i < count; i++)
            {
                term11 = count;
                term12 += xData[i];
                term13 += Math.Pow(xData[i], 2);
                term21 += xData[i];
                term22 += Math.Pow(xData[i], 2);
                term23 += Math.Pow(xData[i], 3);
                term31 += Math.Pow(xData[i], 2);
                term32 += Math.Pow(xData[i], 3);
                term33 += Math.Pow(xData[i], 4);
                termY1 += Math.Log(Math.Abs(yData[i]));
                termY2 += Math.Log(Math.Abs(yData[i])) * xData[i];
                termY3 += Math.Log(Math.Abs(yData[i])) * xData[i] * xData[i];
            }

            double delta = term11 * term22 * term33 + term21 * term32 * term13 + term12 * term23 * term31
                            - term13 * term22 * term31 - term12 * term21 * term33 - term11 * term23 * term32;

            double delta_x1 = termY1 * term22 * term33 + termY2 * term32 * term13 + term12 * term23 * termY3
                            - term13 * term22 * termY3 - term12 * termY2 * term33 - termY1 * term23 * term32;

            double delta_x2 = term11 * termY2 * term33 + term21 * termY3 * term13 + termY1 * term23 * term31
                            - term13 * termY2 * term31 - termY1 * term21 * term33 - term11 * term23 * termY3;

            double delta_x3 = term11 * term22 * termY3 + term21 * term32 * termY1 + term12 * termY2 * term31
                            - termY1 * term22 * term31 - term12 * term21 * termY3 - term11 * termY2 * term32;

            double a = delta_x1 / delta;
            double b = delta_x2 / delta;
            double c = delta_x3 / delta;

            double mu = -1 * b / (2 * c);
            double sigma = Math.Sqrt(-1 / (2 * c));
            double A = Math.Exp(a - b * b / (4 * c));

            gaussianParameter.mu = mu;
            gaussianParameter.sigma = sigma;
            gaussianParameter.amp = A;
        }

        #endregion


        #region Gaussian Fittin  GUOS_ALGORITHM

        void GUOS_ALGORITHM(List<double> xData, List<double> yData,
                                    out GaussianParameter gaussianParameter)
        {
            gaussianParameter = new GaussianParameter();

            double term11 = 0;
            double term12 = 0; 
            double term13 = 0;
            double term21 = 0;
            double term22 = 0;
            double term23 = 0;
            double term31 = 0;
            double term32 = 0;
            double term33 = 0;
            double termY1 = 0;
            double termY2 = 0;
            double termY3 = 0;

            int count = yData.Count;

            for (int i = 0; i < count; i++)
            {
                term11 += Math.Pow(yData[i], 2);
                term12 += xData[i] * Math.Pow(yData[i], 2);
                term13 += Math.Pow(xData[i], 2) * Math.Pow(yData[i], 2);
                term21 += xData[i] * Math.Pow(yData[i], 2);
                term22 += Math.Pow(xData[i], 2) * Math.Pow(yData[i], 2);
                term23 += Math.Pow(xData[i], 3) * Math.Pow(yData[i], 2);
                term31 += Math.Pow(xData[i], 2) * Math.Pow(yData[i], 2);
                term32 += Math.Pow(xData[i], 3) * Math.Pow(yData[i], 2);
                term33 += Math.Pow(xData[i], 4) * Math.Pow(yData[i], 2);
                termY1 += Math.Log(Math.Abs(yData[i])) * Math.Pow(yData[i], 2);
                termY2 += Math.Log(Math.Abs(yData[i])) * xData[i] * Math.Pow(yData[i], 2);
                termY3 += Math.Log(Math.Abs(yData[i])) * xData[i] * xData[i] * Math.Pow(yData[i], 2);
            }

            double delta = term11 * term22 * term33 + term21 * term32 * term13 + term12 * term23 * term31
                            - term13 * term22 * term31 - term12 * term21 * term33 - term11 * term23 * term32;

            double delta_x1 = termY1 * term22 * term33 + termY2 * term32 * term13 + term12 * term23 * termY3
                            - term13 * term22 * termY3 - term12 * termY2 * term33 - termY1 * term23 * term32;

            double delta_x2 = term11 * termY2 * term33 + term21 * termY3 * term13 + termY1 * term23 * term31
                            - term13 * termY2 * term31 - termY1 * term21 * term33 - term11 * term23 * termY3;

            double delta_x3 = term11 * term22 * termY3 + term21 * term32 * termY1 + term12 * termY2 * term31
                            - termY1 * term22 * term31 - term12 * term21 * termY3 - term11 * termY2 * term32;

            double a = delta_x1 / delta;
            double b = delta_x2 / delta;
            double c = delta_x3 / delta;

            double mu = -1 * b / (2 * c);
            double sigma = Math.Sqrt(-1 / (2 * c));
            double A = Math.Exp(a - b * b / (4 * c));

            gaussianParameter.mu = mu;
            gaussianParameter.sigma = sigma;
            gaussianParameter.amp = A;
        }

        #endregion

        #region Gaussian Fitting TwoPoint
        void InputConditionUsingTwoPoint(List<double> xData, List<double> yData, double condition, 
                                    out GaussianParameter gaussianParameter)
        {
            double Max = yData.Max();
            double Target = Max * condition / 100;
            int inderxMax = FindIndex(Max, yData);

            double xLeft = FindLeftPoint(Target, inderxMax, xData, yData);
            double xRight = FindRightPoint(Target, inderxMax, xData, yData);

            double amp = Max;
            double mu = (xRight - xLeft) / 2;
            double a1 = Math.Pow(xRight - mu, 2);
            double a2 = Math.Pow(xLeft - mu, 2);

            double Term1 = -1 * a1 + a2;
            double Term2 = Math.Log(1);
            double sigma = Math.Sqrt(Term1 / (2 * Term2));

            gaussianParameter = new GaussianParameter();
            gaussianParameter.amp = amp;
            gaussianParameter.mu = mu;
            gaussianParameter.sigma = sigma;
        }

        int FindIndex(double Target, List<double> yData)
        {
            int index = 0;

            for (int i = 0; i < yData.Count; i++)
            {
                if (yData[i] == Target)
                {
                    index = i; 
                    break;
                }
            }

            return index;
        }

        double FindLeftPoint(double Target, int inderxMax, List<double> xData, List<double> yData)
        {
            double point = 0;

            for (int n = 0; n < inderxMax - 1; n++)
            {
                int i = inderxMax - n;

                if (yData[i] == Target)
                {
                    point = xData[i];
                    break;
                }
                else if (yData[i] > Target && Target > yData[i - 1])
                {
                    double a = (yData[i] - yData[i - 1]) / (xData[i] - xData[i - 1]);
                    double b = (yData[i - 1] * xData[i] - yData[i] * xData[i - 1]) / (xData[i] - xData[i - 1]);

                    point = (Target - b) / a;
                    break;
                }
            }

            return point;
        }

        double FindRightPoint(double Target, int inderxMax, List<double> xData, List<double> yData)
        {
            double point = 0;

            for (int i = inderxMax; i < yData.Count; i++)
            {
                if (yData[i] == Target)
                {
                    point = xData[i];
                    break;
                }
                else if (yData[i] < Target && Target < yData[i - 1])
                {
                    double a = (yData[i] - yData[i - 1]) / (xData[i] - xData[i - 1]);
                    double b = (yData[i - 1] * xData[i] - yData[i] * xData[i - 1]) / (xData[i] - xData[i - 1]);

                    point = (Target - b) / a;
                    break;
                }
            }

            return point;
        }
        #endregion

        #region Gaussian Function
        void CalculateGaussian(List<double> xData,
                                GaussianParameter gaussianParameter, out List<double> yDataG)
        {
            int dataCount = FFData.xFFV.Count;
            yDataG = new List<double>();

            for (int i = 0; i < dataCount; i++)
            {
                double value = GaussianFunction(xData[i], gaussianParameter);
                yDataG.Add(value);
            }
        }

        double GaussianFunction(double x, GaussianParameter gaussianParameter)
        {
            double mean = gaussianParameter.mu;
            double stdDev = gaussianParameter.sigma;
            double amplitude = gaussianParameter.amp;

            return amplitude * Math.Exp(-Math.Pow((x - mean) / stdDev, 2) / 2);
        }
        #endregion

        #region 數據讀取與輸出
        void LoadData(string dataPath, out FarFieldData FFData)
        {
            FFData = new FarFieldData();

            StreamReader streamReader = new StreamReader(dataPath);
            streamReader.ReadLine(); // Title

            while (streamReader.EndOfStream == false)
            {
                string sLine = streamReader.ReadLine();
                string[] array = sLine.Split(',');

                FFData.xFFH.Add(Convert.ToDouble(array[0]));
                FFData.yFFH.Add(Convert.ToDouble(array[1]));
                FFData.xFFV.Add(Convert.ToDouble(array[2]));
                FFData.yFFV.Add(Convert.ToDouble(array[3]));
            }

            streamReader.Close();
        }

        void SaveData(string dataPath, FarFieldData FFData, List<double> yDataG)
        {
            int count  = yDataG.Count;

            StreamWriter streamWriter = new StreamWriter(dataPath);
            streamWriter.WriteLine("x[D],y[], y_G[]");

            for (int i = 0; i < count; i++)
            {
                streamWriter.WriteLine($"{FFData.xFFV[i]},{FFData.yFFV[i]},{yDataG[i]}");
            }

            streamWriter.Close();
        }

        void SaveData(string dataPath, FarFieldData FFData, List<double> yHDataG, List<double> yVDataG)
        {
            int count = FFData.xFFH.Count;

            StreamWriter streamWriter = new StreamWriter(dataPath);
            streamWriter.WriteLine("x[D],yH[], yHG[],x[D],yV[], yVG[]");

            for (int i = 0; i < count; i++)
            {
                streamWriter.WriteLine($"{FFData.xFFH[i]},{FFData.yFFH[i]},{yHDataG[i]}, {FFData.xFFV[i]},{FFData.yFFV[i]},{yVDataG[i]}");
            }

            streamWriter.Close();
        }

        void SaveData(string dataPath, FarFieldData FFData)
        {
            StreamWriter streamWriter = new StreamWriter(dataPath);
            streamWriter.WriteLine("x_h[D],y_h[], x_v[D],y_v[]");

            double xValue = -45;

            for (int i = 0; i <901; i++)
            {
                if (FFData.xFFV[i] == xValue)
                {
                    streamWriter.WriteLine($"{xValue},{FFData.yFFH[i]},{xValue},{FFData.yFFV[i]}");
                    xValue++;
                }
                
            }

            streamWriter.Close();
        }
        #endregion
    }

    class GaussianParameter
    {
        public double amp = 0;
        public double mu = 0;
        public double sigma = 0;
    }

    class FarFieldData
    {
        public List<double> xFFH = new List<double>();
        public List<double> yFFH = new List<double>();
        public List<double> xFFV = new List<double>();
        public List<double> yFFV = new List<double>();
    }
}
