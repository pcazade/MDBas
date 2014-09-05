#include "mainwindow.h"
#include "crystal_to_box.h"

MainWindow::MainWindow()
{
    // parameters input
    a = new QLineEdit("0.0");
    b = new QLineEdit("0.0");
    c = new QLineEdit("0.0");
    alpha = new QLineEdit("0.0");
    beta  = new QLineEdit("0.0");
    gamma = new QLineEdit("0.0");

    QLabel *name_a = new QLabel("Length a");
    QLabel *name_b = new QLabel("Length b");
    QLabel *name_c = new QLabel("Length c");
    QLabel *name_alph = new QLabel("Angle alpha");
    QLabel *name_beta = new QLabel("Angle beta");
    QLabel *name_gamm = new QLabel("Angle gamma");

    // organisation with layout
    QGridLayout *latticeParamsLayout = new QGridLayout;

    latticeParamsLayout->addWidget(name_a,0,0);
    latticeParamsLayout->addWidget(name_b,0,1);
    latticeParamsLayout->addWidget(name_c,0,2);
    latticeParamsLayout->addWidget(a,1,0);
    latticeParamsLayout->addWidget(b,1,1);
    latticeParamsLayout->addWidget(c,1,2);

    latticeParamsLayout->addWidget(name_alph,2,0);
    latticeParamsLayout->addWidget(name_beta,2,1);
    latticeParamsLayout->addWidget(name_gamm,2,2);
    latticeParamsLayout->addWidget(alpha,3,0);
    latticeParamsLayout->addWidget(beta,3,1);
    latticeParamsLayout->addWidget(gamma,3,2);

    // a distinct box for storing previous things
    QGroupBox *latticeParamsGroup = new QGroupBox("Lattice parameters: ");
    latticeParamsGroup->setLayout(latticeParamsLayout);

    //----------------------------------------------------------------------
    //output zone : matrix
    m1 = new QLineEdit("0.0");
    m2 = new QLineEdit("0.0");
    m3 = new QLineEdit("0.0");
    m4 = new QLineEdit("0.0");
    m5 = new QLineEdit("0.0");
    m6 = new QLineEdit("0.0");
    m7 = new QLineEdit("0.0");
    m8 = new QLineEdit("0.0");
    m9 = new QLineEdit("0.0");

    // organisation with a grid layout
    QGridLayout *matrixLayout = new QGridLayout;
    matrixLayout->addWidget(m1,0,0);
    matrixLayout->addWidget(m2,0,1);
    matrixLayout->addWidget(m3,0,2);
    matrixLayout->addWidget(m4,1,0);
    matrixLayout->addWidget(m5,1,1);
    matrixLayout->addWidget(m6,1,2);
    matrixLayout->addWidget(m7,2,0);
    matrixLayout->addWidget(m8,2,1);
    matrixLayout->addWidget(m9,2,2);

    // a distinct box for this Matrix
    QGroupBox *matrixGroup = new QGroupBox("Output Matrix ");
    matrixGroup->setLayout(matrixLayout);

    //----------------------------------------------------------------------
    //buttons
    launchButton = new QPushButton("Obtain matrix",this);
    resetButton = new QPushButton("Reset",this);
    exitButton = new QPushButton("Exit",this);

    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addWidget(launchButton);
    buttonsLayout->addWidget(resetButton);
    buttonsLayout->addWidget(exitButton);

    //----------------------------------------------------------------------

    // general layout of window : we add here the previous group boxes, etc ...
    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(latticeParamsGroup);
    mainLayout->addLayout(buttonsLayout);
    mainLayout->addWidget(matrixGroup);

    // general settings for this main window
    this->setLayout(mainLayout);
    this->setWindowTitle("Lattice param to Matrix");
    this->setFixedSize(400,500);

    // connect buttons' signals to something else ...
    connect(exitButton,SIGNAL(clicked()),qApp,SLOT(quit()));
    connect(launchButton,SIGNAL(clicked()),this,SLOT(computeMatrix()));
    connect(resetButton,SIGNAL(clicked()),this,SLOT(resetAll()));
}


void MainWindow::computeMatrix()
{
    real lattice[6]= {0.0} , matrix[9]= {0.0};
    int errorCode=0;

    lattice[0] = this->a->text().toDouble();
    lattice[1] = this->b->text().toDouble();
    lattice[2] = this->c->text().toDouble();
    lattice[3] = this->alpha->text().toDouble();
    lattice[4] = this->beta->text().toDouble();
    lattice[5] = this->gamma->text().toDouble();

    try
    {
        errorCode = lattice_to_cryst(lattice,matrix);

        if (errorCode!=0)
            throw errorCode;
    }
    catch(int e)
    {
        QString lapack_error = "The subroutine DSYEV of LAPACK returned the following error code : " + QString::number(e) + "\n" + "Please check LAPACK documentation.";
        QMessageBox::critical(this, "LAPACK ERROR !",lapack_error);
        return;
    }

    this->m1->setText(QString::number(matrix[0]));
    this->m2->setText(QString::number(matrix[1]));
    this->m3->setText(QString::number(matrix[2]));

    this->m4->setText(QString::number(matrix[3]));
    this->m5->setText(QString::number(matrix[4]));
    this->m6->setText(QString::number(matrix[5]));

    this->m7->setText(QString::number(matrix[6]));
    this->m8->setText(QString::number(matrix[7]));
    this->m9->setText(QString::number(matrix[8]));

}

void MainWindow::resetAll()
{
    this->a->setText(QString::number(0.0));
    this->b->setText(QString::number(0.0));
    this->c->setText(QString::number(0.0));

    this->alpha->setText(QString::number(0.0));
    this->beta->setText(QString::number(0.0));
    this->gamma->setText(QString::number(0.0));

    this->m1->setText(QString::number(0.0));
    this->m2->setText(QString::number(0.0));
    this->m3->setText(QString::number(0.0));

    this->m4->setText(QString::number(0.0));
    this->m5->setText(QString::number(0.0));
    this->m6->setText(QString::number(0.0));

    this->m7->setText(QString::number(0.0));
    this->m8->setText(QString::number(0.0));
    this->m9->setText(QString::number(0.0));
}
