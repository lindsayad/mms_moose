#include "ArticunoApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

#include "CoupledGradientBC.h"

template<>
InputParameters validParams<ArticunoApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

ArticunoApp::ArticunoApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  ArticunoApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  ArticunoApp::associateSyntax(_syntax, _action_factory);
}

ArticunoApp::~ArticunoApp()
{
}

// External entry point for dynamic application loading
extern "C" void ArticunoApp__registerApps() { ArticunoApp::registerApps(); }
void
ArticunoApp::registerApps()
{
  registerApp(ArticunoApp);
}

// External entry point for dynamic object registration
extern "C" void ArticunoApp__registerObjects(Factory & factory) { ArticunoApp::registerObjects(factory); }
void
ArticunoApp::registerObjects(Factory & factory)
{
  registerBoundaryCondition(CoupledGradientBC);
}

// External entry point for dynamic syntax association
extern "C" void ArticunoApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { ArticunoApp::associateSyntax(syntax, action_factory); }
void
ArticunoApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
