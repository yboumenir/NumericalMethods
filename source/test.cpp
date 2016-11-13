#include <QApplication>
#include "../include/numericalmatrix.h"
#include "../include/numericalmethods.h"



#include "qwt_plot.h"
#include "qwt_plot_curve.h"
#include "qwt_plot_grid.h"
#include "qwt_symbol.h"
#include "qwt_legend.h"
#include <iostream>

struct plotStruct{


public:

    plotStruct(QString plotName = "Convergance Plot"){

        plot = new QwtPlot();
        grid = new QwtPlotGrid();
        curve = new QwtPlotCurve();
        symbol = new QwtSymbol(QwtSymbol::Ellipse, QBrush( Qt::green ), QPen( Qt::green, 1 ), QSize( 4, 4 ));
        points = new QPolygonF();


        plot->setTitle( plotName );
        plot->setCanvasBackground( Qt::white );
        plot->setAxisScale( QwtPlot::yLeft, 0.0, 10.0 );
        plot->insertLegend( new QwtLegend() );


        grid->attach( plot );

        curve->setTitle( "Solutions" );


        plot->resize( 600, 400 );



    }


    void addPoints(QPolygonF points, QString plotTitle = "plot_"){
        static int curve_id;


        QwtPlotCurve *curve = new QwtPlotCurve;


        if(!plotTitle.compare("plot_",Qt::CaseInsensitive))
        {
            curve_id++;
            std::string tmp ("plot_");
            tmp += std::to_string(curve_id);
            curve->setTitle(QString::fromStdString(tmp));
        }
        else
        {
            curve->setTitle(plotTitle);
        }

        int penid = Qt::color0 + curve_id;

        std::cout <<penid;

        curve->setPen( penid, 2 ), curve->setRenderHint( QwtPlotItem::RenderAntialiased, true );

        curve->setSymbol( symbol );

        curve->setSamples( points );
        curve->attach( plot );

    }

    void show(){
        plot->show();
    }



private:
    QwtPlot *plot;
    QwtPlotGrid *grid;
    QwtPlotCurve *curve;
    QwtSymbol *symbol;
    QPolygonF *points;


};


int xmain( int argc, char **argv )
{


    plotStruct plt;



    QPolygonF point;
    for(double i=0.0;i<6.28;i = i + 0.0628 ){
        point<<QPointF(i,std::sin(i)+1);

    }
    plt.addPoints(point);


    plt.show();

}

struct Functions{

    plotStruct plt;

    static double foo(double x){
        return pow(cos(x),2)+sin(x);
    }

    // need to define the set of functions
    NumericalMatrix<double> F;
    // next need to define the jacobians of F
    NumericalMatrix<double> JJ;
    // internal TODO: needs to be internal of solver
    NumericalMatrix<double> iJJ;
    // The guess vector
    NumericalMatrix<double> X;
    void initializeFunctions(NumericalMatrix<double> & X,
                             NumericalMatrix<double>& JJ,
                             NumericalMatrix<double> & F,
                             NumericalMatrix<double> & iJJ){
        X.createVector({2,1});
        JJ.createMatrix(2,2);
        F.createVector({0,0});
        iJJ.createMatrix(2,2);
    }

    // error
    double eps = .000001;

    void initializeJacobian(NumericalMatrix<double> & JJ){
        // jacobian with initial guess
        JJ.setValue(0,0,2*X.getElement(0));
        JJ.setValue(0,1,2*X.getElement(1));
        JJ.setValue(1,0,2*X.getElement(1));
        JJ.setValue(1,1,2*X.getElement(0));
    }
    // the sys of non eqs.
    //f1 = (x1)^2 +(x2)^2 - 50 =0
    //f2 = (x1)*(x2) -25=0
    void computeFunction(NumericalMatrix<double> & F){
        this->F.setValue(0,pow(this->X.getElement(0), 2) + pow(this->X.getElement(1), 2) - 50);
        this->F.setValue(1,this->X.getElement(0) * this->X.getElement(1) - 25);
    }

    void computeInverseJacobian(NumericalMatrix<double> JJ, NumericalMatrix<double> & iJJ){
        iJJ = NumericalMethods<double>::inverse(JJ);

    }

    bool computeNewGuess(NumericalMatrix<double> &X, int i){
        X.setValue(0,X.getElement(0) - (this->iJJ[0][0] * this->F.getElement(0) + this->iJJ[0][1] * this->F.getElement(1)));
        X.setValue(1,X.getElement(1) - (this->iJJ[1][0] * this->F.getElement(1) + this->iJJ[1][1] * this->F.getElement(1)));

        static QPolygonF point1, point2;

        point1 << QPointF(i,X.getElement(0));
        point2 << QPointF(i,X.getElement(1));
        std::cout << "x(1) = " << X.getElement(0) << ", x(2) = " << X.getElement(1) << std::endl;
        if ((abs(pow(X.getElement(0), 2) + pow(X.getElement(1), 2) - 50) < eps) && abs(X.getElement(0) * X.getElement(1) - 25)<eps)
        {
            printf("Solution converged at N = %3.3i to X[0] = %3.3f and X[1] = %3.3f\r\n",i, X.getElement(0), X.getElement(1));
            plt.addPoints(point1);
            plt.addPoints(point2);
            plt.show();
            return 1;
        }

    }




};

#include <functional>
#include <memory>

int main(int argc, char **argv)
{
    QApplication a(argc,argv);
    Functions funcs;

    try {
        NumericalMatrix<double>  mat;
        int row = 3;
        int col = 3;

        mat.createMatrix(row,col);

        mat.pushColumns(0,{7,2,1});
        mat.pushColumns(1,{0,3,-1});
        mat.pushColumns(2,{-3,4,-2});

        std::cout << "A: " << std::endl;
        mat.printMatrix();

        NumericalMethods<double>::transpose(mat);

        std::cout << "det: " << NumericalMethods<double>::determinant(mat.dimension(),mat) << std::endl;

        std::cout << "inv: "<< std::endl;
        NumericalMethods<double>::inverse(mat).printMatrix();

        std::cout << "newton solver" << std::endl;
        NumericalMethods<double>::newtonsMethod(funcs);

        std::cout << "integrating cos(x) from 0 to 3.14/2  :" << std::endl;
        std::cout << NumericalMethods<double>::trap(Functions::foo,0,M_PI/2);
        std::cout << std::endl;

//        intStack.createVector({1,2});

//        // manipulate int stack
//        intStack.push(7);
//        intStack.printVector();
//        std::cout << intStack.top() <<std::endl;


    }
    catch (std::exception const& ex) {
        std::cerr << "Exception: " << ex.what() <<std::endl;
        return -1;
    }

    return a.exec();
}

