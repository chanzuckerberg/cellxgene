// src/components/AuthN.js

import React from "react";
import { useAuth0 } from "../../react-auth0-spa";

const AuthN = () => {
  const { isAuthenticated, loginWithRedirect, logout } = useAuth0();

  return (
    <div>
      {!isAuthenticated && (
        <button type="button" onClick={() => loginWithRedirect({})}>
          Log in
        </button>
      )}

      {isAuthenticated && (
        <button type="button" onClick={() => logout()}>
          Log out
        </button>
      )}
    </div>
  );
};

export default AuthN;
