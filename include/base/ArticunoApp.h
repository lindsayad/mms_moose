#ifndef ARTICUNOAPP_H
#define ARTICUNOAPP_H

#include "MooseApp.h"

class ArticunoApp;

template<>
InputParameters validParams<ArticunoApp>();

class ArticunoApp : public MooseApp
{
public:
  ArticunoApp(InputParameters parameters);
  virtual ~ArticunoApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* ARTICUNOAPP_H */
