import { Position, Toaster, Intent } from "@blueprintjs/core";

/** Singleton toaster instance. Create separate instances for different options. */

const ErrorToastTopCenter = Toaster.create({
  className: "recipe-toaster",
  position: Position.TOP
});

/*
A "user" error - eg, bad input
*/
export const postUserErrorToast = message =>
  ErrorToastTopCenter.show({ message, intent: Intent.WARNING });

/*
A toast the user must dismiss manually, because they need to act on its information,
ie., 8 bulk add genes out of 40 were bad. Manually see which ones and fix.
*/
export const keepAroundErrorToast = message =>
  ErrorToastTopCenter.show({ message, timeout: 0, intent: Intent.WARNING });

/*
a hard network error
*/
export const postNetworkErrorToast = message =>
  ErrorToastTopCenter.show({
    message,
    timeout: 30000,
    intent: Intent.DANGER
  });
