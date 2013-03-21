#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui>

class MainWindow : public QWidget
{
    Q_OBJECT

public:
    //constructor
    MainWindow();

private slots:
    void computeMatrix();
    void resetAll();

private:
    // for lattice params
    QLineEdit *a , *b , *c;
    QLineEdit *alpha , *beta, *gamma;

    // for matrix ; m1 to m3 for first row, etc...
    QLineEdit *m1,*m2,*m3,*m4,*m5;
    QLineEdit *m6,*m7,*m8,*m9;

    // button for launching calculation
    QPushButton *launchButton;
    QPushButton *resetButton;
    QPushButton *exitButton;

};

#endif // MAINWINDOW_H
