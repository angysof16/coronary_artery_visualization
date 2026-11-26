#include <vtkSmartPointer.h>
#include <vtkAlgorithmOutput.h>
#include <vtkMetaImageReader.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkExtractVOI.h>
#include <vtkMarchingCubes.h>
#include <vtkPointLocator.h>
#include <vtkClipPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkTubeFilter.h>
#include <vtkColorTransferFunction.h>
#include <vtkMath.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cout << "Uso: " << argv[0]
              << " image.mhd linea_central.vtk [umbral] [dist_mm] [umbral_corazon] [umbral_calcio]\n"
              << "  umbral         : valor de la isosuperficie. Si se omite o es 0,\n"
              << "                   se usa la mitad de la intensidad en la arteria.\n"
              << "  dist_mm        : si es mayor que 0, aisla la arteria, calcula el radio\n"
              << "                   y colorea el tubo por estrechamiento.\n"
              << "                   Si se omite o es 0, muestra el checkpoint 0.\n"
              << "  umbral_corazon : si es mayor que 0, dibuja ademas el corazon\n"
              << "                   (sangre con contraste) con ese umbral.\n"
              << "                   Si se omite o es 0, solo se ve la arteria.\n"
              << "  umbral_calcio  : si es mayor que 0 (y dist_mm > 0), dibuja las\n"
              << "                   calcificaciones cerca de la linea central.\n"
              << "                   Debe ser claramente mayor que la intensidad de la arteria.\n"
              << "Ejemplos:\n"
              << "  " << argv[0] << " image.mhd Team_LKEB_ManuallyCorrected_vessel00.vtk\n"
              << "  " << argv[0] << " image.mhd Team_LKEB_ManuallyCorrected_vessel00.vtk 0 6\n"
              << "  " << argv[0] << " image.mhd Team_LKEB_ManuallyCorrected_vessel00.vtk 0 6 250\n"
              << "  " << argv[0] << " image.mhd Team_LKEB_ManuallyCorrected_vessel00.vtk 0 6 0 600\n";
    return 1;
  }

  // ------------------------------------------------------------------
  // Lectura del volumen y de la linea central
  // ------------------------------------------------------------------
  auto reader = vtkSmartPointer<vtkMetaImageReader>::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  vtkImageData* img = reader->GetOutput();

  auto clReader = vtkSmartPointer<vtkPolyDataReader>::New();
  clReader->SetFileName(argv[2]);
  clReader->Update();
  vtkPolyData* cl = clReader->GetOutput();

  // Comparar sistemas de coordenadas
  double ib[6], cb[6];
  img->GetBounds(ib);
  cl->GetBounds(cb);
  std::cout << "Volumen bounds: " << ib[0] << " " << ib[1] << " " << ib[2] << " "
            << ib[3] << " " << ib[4] << " " << ib[5] << "\n";
  std::cout << "Linea bounds:   " << cb[0] << " " << cb[1] << " " << cb[2] << " "
            << cb[3] << " " << cb[4] << " " << cb[5] << "\n";

  // Intensidad del volumen en el punto medio de la linea central
  double pMid[3], sp[3], org[3];
  int ext[6];
  cl->GetPoint(cl->GetNumberOfPoints() / 2, pMid);
  img->GetSpacing(sp);
  img->GetOrigin(org);
  img->GetExtent(ext);
  int i = std::lround((pMid[0] - org[0]) / sp[0]);
  int j = std::lround((pMid[1] - org[1]) / sp[1]);
  int k = std::lround((pMid[2] - org[2]) / sp[2]);
  double lumen = img->GetScalarComponentAsDouble(i, j, k, 0);
  std::cout << "Intensidad en la arteria: " << lumen << "\n";

  double umbral = (argc > 3 && std::atof(argv[3]) > 0) ? std::atof(argv[3]) : lumen * 0.5;
  double dist = (argc > 4) ? std::atof(argv[4]) : 0.0;
  // El corazon solo se dibuja si se da un umbral para el en el quinto argumento
  double umbralCorazon = (argc > 5) ? std::atof(argv[5]) : 0.0;
  bool conCorazon = (dist > 0 && umbralCorazon > 0);
  // Las calcificaciones solo se dibujan si se da un umbral en el sexto argumento
  double umbralCalcio = (argc > 6) ? std::atof(argv[6]) : 0.0;
  bool conCalcio = (dist > 0 && umbralCalcio > 0);
  std::cout << "Umbral: " << umbral << "\n";
  if (conCalcio && umbralCalcio <= lumen)
  {
    std::cout << "Aviso: el umbral de calcio (" << umbralCalcio
              << ") no supera la intensidad de la arteria (" << lumen
              << "). Probablemente se vera toda la arteria como calcio.\n";
  }

  // Pasa una caja en mm (la de la linea central mas un margen) a indices de voxel
  auto cajaVoxel = [&](double margen, int voi[6])
  {
    for (int d = 0; d < 3; ++d)
    {
      int lo = (int)std::floor((cb[2 * d] - margen - org[d]) / sp[d]);
      int hi = (int)std::ceil((cb[2 * d + 1] + margen - org[d]) / sp[d]);
      voi[2 * d] = std::max(lo, ext[2 * d]);
      voi[2 * d + 1] = std::min(hi, ext[2 * d + 1]);
    }
  };

  // ------------------------------------------------------------------
  // Marching Cubes. Con dist > 0 se aplica antes sobre un volumen recortado.
  // ------------------------------------------------------------------
  auto mc = vtkSmartPointer<vtkMarchingCubes>::New();
  mc->SetValue(0, umbral);
  mc->ComputeNormalsOn();

  auto crop = vtkSmartPointer<vtkExtractVOI>::New();
  if (dist > 0)
  {
    int voi[6];
    cajaVoxel(15.0, voi);
    crop->SetInputData(img);
    crop->SetVOI(voi);
    mc->SetInputConnection(crop->GetOutputPort());
  }
  else
  {
    mc->SetInputData(img);
  }

  vtkAlgorithmOutput* salida = mc->GetOutputPort();

  // Objetos que deben seguir vivos hasta el render
  auto locator = vtkSmartPointer<vtkPointLocator>::New();
  auto locSup = vtkSmartPointer<vtkPointLocator>::New();
  auto conDist = vtkSmartPointer<vtkPolyData>::New();
  auto clip = vtkSmartPointer<vtkClipPolyData>::New();
  auto conn = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  auto sArr = vtkSmartPointer<vtkFloatArray>::New();
  sArr->SetName("estrechamiento");
  std::vector<double> radio;

  // ------------------------------------------------------------------
  // Paso 1 y B. Aislar la arteria y calcular el radio por punto
  // ------------------------------------------------------------------
  if (dist > 0)
  {
    mc->Update();
    vtkPolyData* cruda = mc->GetOutput();
    std::cout << "Triangulos de la superficie cruda: " << cruda->GetNumberOfPolys() << "\n";

    // Distancia de cada vertice al punto mas cercano de la linea central
    locator->SetDataSet(cl);
    locator->BuildLocator();
    auto d = vtkSmartPointer<vtkFloatArray>::New();
    d->SetName("distancia");
    for (vtkIdType n = 0; n < cruda->GetNumberOfPoints(); ++n)
    {
      double p[3], q[3];
      cruda->GetPoint(n, p);
      cl->GetPoint(locator->FindClosestPoint(p), q);
      d->InsertNextValue((float)std::sqrt(vtkMath::Distance2BetweenPoints(p, q)));
    }
    conDist->ShallowCopy(cruda);
    conDist->GetPointData()->SetScalars(d);

    // Se conserva lo que esta a menos de dist mm
    clip->SetInputData(conDist);
    clip->SetValue(dist);
    clip->InsideOutOn();
    clip->GenerateClipScalarsOff();
    clip->Update();
    std::cout << "Triangulos tras el recorte: " << clip->GetOutput()->GetNumberOfPolys() << "\n";

    // Se deja la pieza conectada mas cercana al punto medio de la linea
    conn->SetInputConnection(clip->GetOutputPort());
    conn->SetExtractionModeToClosestPointRegion();
    conn->SetClosestPoint(pMid);
    conn->Update();
    std::cout << "Triangulos tras la conectividad: " << conn->GetOutput()->GetNumberOfPolys() << "\n";

    // Radio en cada punto de la linea central, distancia al vertice mas cercano de la pared
    vtkPolyData* arteria = conn->GetOutput();
    locSup->SetDataSet(arteria);
    locSup->BuildLocator();
    for (vtkIdType n = 0; n < cl->GetNumberOfPoints(); ++n)
    {
      double p[3], q[3];
      cl->GetPoint(n, p);
      arteria->GetPoint(locSup->FindClosestPoint(p), q);
      radio.push_back(std::sqrt(vtkMath::Distance2BetweenPoints(p, q)));
    }

    // Suavizado, cada radio pasa a ser el promedio de sus 5 vecinos
    std::vector<double> suave(radio.size());
    for (int n = 0; n < (int)radio.size(); ++n)
    {
      double suma = 0;
      int cuenta = 0;
      for (int kk = -2; kk <= 2; ++kk)
      {
        int m = n + kk;
        if (m < 0 || m >= (int)radio.size()) continue;
        suma += radio[m];
        cuenta++;
      }
      suave[n] = suma / cuenta;
    }
    radio = suave;

    // Radio de referencia, el percentil 90 de todos los radios
    std::vector<double> ordenado = radio;
    std::sort(ordenado.begin(), ordenado.end());
    double rRef = ordenado[(size_t)(0.9 * (ordenado.size() - 1))];

    // Estrechamiento en cada punto, entre 0 y 1
    double sMax = 0;
    size_t iMax = 0;
    for (size_t n = 0; n < radio.size(); ++n)
    {
      double s = std::min(1.0, std::max(0.0, 1.0 - radio[n] / rRef));
      sArr->InsertNextValue((float)s);
      if (s > sMax) { sMax = s; iMax = n; }
    }
    std::cout << "Radio min " << ordenado.front() << " mm, max " << ordenado.back()
              << " mm, referencia (p90) " << rRef << " mm\n"
              << "Mayor estrechamiento " << sMax << " en el punto " << iMax
              << " de " << radio.size() << "\n";

    salida = conn->GetOutputPort();
  }

  // ------------------------------------------------------------------
  // Calcificaciones (opcional). Isosuperficie alta cerca de la linea central.
  // Reusa el locator de la linea central construido en el bloque anterior.
  // ------------------------------------------------------------------
  auto cropC = vtkSmartPointer<vtkExtractVOI>::New();
  auto mcC = vtkSmartPointer<vtkMarchingCubes>::New();
  auto conDistC = vtkSmartPointer<vtkPolyData>::New();
  auto clipC = vtkSmartPointer<vtkClipPolyData>::New();

  if (conCalcio)
  {
    int voiC[6];
    cajaVoxel(15.0, voiC);
    cropC->SetInputData(img);
    cropC->SetVOI(voiC);

    mcC->SetInputConnection(cropC->GetOutputPort());
    mcC->SetValue(0, umbralCalcio);
    mcC->ComputeNormalsOn();
    mcC->Update();
    vtkPolyData* calcioCrudo = mcC->GetOutput();

    // Distancia de cada vertice a la linea central
    auto dC = vtkSmartPointer<vtkFloatArray>::New();
    dC->SetName("distancia");
    for (vtkIdType n = 0; n < calcioCrudo->GetNumberOfPoints(); ++n)
    {
      double p[3], q[3];
      calcioCrudo->GetPoint(n, p);
      cl->GetPoint(locator->FindClosestPoint(p), q);
      dC->InsertNextValue((float)std::sqrt(vtkMath::Distance2BetweenPoints(p, q)));
    }
    conDistC->ShallowCopy(calcioCrudo);
    conDistC->GetPointData()->SetScalars(dC);

    // Se conserva solo lo que esta a menos de dist mm (descarta huesos, valvulas, etc.)
    clipC->SetInputData(conDistC);
    clipC->SetValue(dist);
    clipC->InsideOutOn();
    clipC->GenerateClipScalarsOff();
    clipC->Update();
    std::cout << "Umbral del calcio: " << umbralCalcio << "\n";
    std::cout << "Triangulos de calcio: " << clipC->GetOutput()->GetNumberOfPolys() << "\n";
  }

  // ------------------------------------------------------------------
  // Corazon de contexto (opcional). Mismo recorrido que la arteria, pero con
  // una caja de 60 mm de margen y sin recorte por distancia.
  // Solo se ve la sangre con contraste (aorta y cavidades conectadas).
  // ------------------------------------------------------------------
  auto cropH = vtkSmartPointer<vtkExtractVOI>::New();
  auto mcH = vtkSmartPointer<vtkMarchingCubes>::New();
  auto connH = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  double cajaVista[6] = {cb[0], cb[1], cb[2], cb[3], cb[4], cb[5]};

  if (conCorazon)
  {
    int voiH[6];
    cajaVoxel(60.0, voiH);
    cropH->SetInputData(img);
    cropH->SetVOI(voiH);

    mcH->SetInputConnection(cropH->GetOutputPort());
    mcH->SetValue(0, umbralCorazon);
    mcH->ComputeNormalsOn();

    connH->SetInputConnection(mcH->GetOutputPort());
    connH->SetExtractionModeToClosestPointRegion();
    connH->SetClosestPoint(pMid);
    connH->Update();
    std::cout << "Umbral del corazon: " << umbralCorazon << "\n";
    std::cout << "Triangulos del corazon: " << connH->GetOutput()->GetNumberOfPolys() << "\n";

    // La camara mira al centro del corazon, no al de la arteria
    connH->GetOutput()->GetBounds(cajaVista);
  }

  // ------------------------------------------------------------------
  // Actores
  // ------------------------------------------------------------------
  auto sMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  sMapper->SetInputConnection(salida);
  sMapper->ScalarVisibilityOff();
  auto sActor = vtkSmartPointer<vtkActor>::New();
  sActor->SetMapper(sMapper);
  sActor->GetProperty()->SetColor(0.85, 0.35, 0.30);
  sActor->GetProperty()->SetOpacity(dist > 0 ? 0.25 : 0.3);

  // Corazon, superficie muy transparente de otro color para no confundirla
  auto hMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  hMapper->SetInputConnection(connH->GetOutputPort());
  hMapper->ScalarVisibilityOff();
  auto hActor = vtkSmartPointer<vtkActor>::New();
  hActor->SetMapper(hMapper);
  hActor->GetProperty()->SetColor(0.80, 0.80, 0.85);
  hActor->GetProperty()->SetOpacity(0.08);

  // Calcificaciones, opacas y de color claro para que contrasten con la pared
  auto calMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  calMapper->SetInputConnection(clipC->GetOutputPort());
  calMapper->ScalarVisibilityOff();
  auto calActor = vtkSmartPointer<vtkActor>::New();
  calActor->SetMapper(calMapper);
  calActor->GetProperty()->SetColor(1.0, 0.95, 0.85);
  calActor->GetProperty()->SetOpacity(1.0);

  // Paso C. Tubo coloreado por estrechamiento (solo si se aislo la arteria)
  auto clS = vtkSmartPointer<vtkPolyData>::New();
  auto tubo = vtkSmartPointer<vtkTubeFilter>::New();
  auto color = vtkSmartPointer<vtkColorTransferFunction>::New();
  auto tMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  auto tActor = vtkSmartPointer<vtkActor>::New();
  auto cMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  auto cActor = vtkSmartPointer<vtkActor>::New();

  // Actor de la linea central. Es el tubo coloreado o la linea amarilla
  vtkActor* lineaActor = nullptr;

  if (dist > 0)
  {
    // Copia de la linea central que lleva el estrechamiento como escalar
    clS->ShallowCopy(cl);
    clS->GetPointData()->SetScalars(sArr);

    // Tubo delgado sobre la linea. El radio es solo visual
    tubo->SetInputData(clS);
    tubo->SetRadius(0.3);
    tubo->SetNumberOfSides(12);

    // Mapa de color viridis, de 0 (sano) a 1 (muy estrecho)
    color->AddRGBPoint(0.00, 0.267, 0.005, 0.329);
    color->AddRGBPoint(0.25, 0.230, 0.322, 0.546);
    color->AddRGBPoint(0.50, 0.128, 0.567, 0.551);
    color->AddRGBPoint(0.75, 0.369, 0.789, 0.383);
    color->AddRGBPoint(1.00, 0.993, 0.906, 0.144);

    tMapper->SetInputConnection(tubo->GetOutputPort());
    tMapper->SetLookupTable(color);
    tMapper->SetScalarRange(0.0, 1.0);
    tMapper->ScalarVisibilityOn();
    tActor->SetMapper(tMapper);
    lineaActor = tActor;
  }
  else
  {
    // Checkpoint 0. Linea central amarilla
    cMapper->SetInputData(cl);
    cActor->SetMapper(cMapper);
    cActor->GetProperty()->SetColor(1.0, 0.9, 0.1);
    cActor->GetProperty()->SetLineWidth(3.0);
    lineaActor = cActor;
  }

  // ------------------------------------------------------------------
  // Una sola vista, con camara oblicua que luego se mueve con el mouse
  // ------------------------------------------------------------------
  auto ren = vtkSmartPointer<vtkRenderer>::New();
  ren->SetBackground(0.25, 0.25, 0.30);
  if (conCorazon) ren->AddActor(hActor);   // primero el corazon
  ren->AddActor(sActor);
  if (conCalcio) ren->AddActor(calActor);  // calcio antes de la linea central
  ren->AddActor(lineaActor);

  double centro[3] = {
    0.5 * (cajaVista[0] + cajaVista[1]),
    0.5 * (cajaVista[2] + cajaVista[3]),
    0.5 * (cajaVista[4] + cajaVista[5])
  };
  vtkCamera* cam = ren->GetActiveCamera();
  cam->SetFocalPoint(centro);
  cam->SetPosition(centro[0] + 1, centro[1] - 1, centro[2] + 1);
  cam->SetViewUp(0, 0, 1);
  // ResetCamera aleja la camara lo necesario para que todo quepa
  ren->ResetCamera();

  auto win = vtkSmartPointer<vtkRenderWindow>::New();
  win->AddRenderer(ren);
  win->SetSize(1000, 800);

  auto iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(win);
  iren->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
  win->Render();
  iren->Start();
}
